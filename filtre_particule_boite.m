% =========================================================================
% FONCTIONS LOCALES (BOX PARTICLE FILTER)
% =========================================================================

function [tracks, detections_associees, next_display_id, next_id] = filtre_particule_boite(tracks, matrice_centres, dt, N_boites, V_modele, erreur_lidar_max, R_cov, next_display_id, next_id)
    
    nb_objets = size(matrice_centres, 1);
    detections_associees = false(1, nb_objets);

    A_int = intval([eye(3), eye(3)*dt, eye(3)*(0.5*dt^2);
                    zeros(3,3), eye(3), eye(3)*dt;
                    zeros(3,3), zeros(3,3), eye(3)]);
    C_int = intval([eye(3), zeros(3,6)]);

    % 1. PRÉDICTION DES BOÎTES
    rayon_max = [2; 2; 2; 1.5; 1.5; 1.5; 1.5; 1.5; 1.5]; % Rayon max autorisé (position, vitesse, accélération)
    for t = 1:numel(tracks)
        % SAUVEGARDE DE LA POSITION GLOBALE (MÉTHODE ENSEMBLISTE)
        tracks(t).X_prec_estime = tracks(t).X_k_estime; 
        
        for b = 1:N_boites
            tracks(t).X_boxes(:, b) = A_int * tracks(t).X_boxes(:, b) + V_modele;
            
            % SATURATION DES BOÎTES (évite l'explosion par wrapping effect)
            centre_b = mid(tracks(t).X_boxes(:, b));
            rayon_b  = min(rad(tracks(t).X_boxes(:, b)), rayon_max);
            tracks(t).X_boxes(:, b) = midrad(centre_b, rayon_b);
        end
        
        centers = zeros(3, N_boites);
        for b = 1:N_boites
            centers(:, b) = mid(tracks(t).X_boxes(1:3, b));
        end
        tracks(t).X_k_estime(1:3) = sum(centers .* repmat(tracks(t).poids, 3, 1), 2);
    end

    % 2. ASSOCIATION ET CORRECTION
    for t = 1:numel(tracks)
        if nb_objets == 0, break; end
        
        centre_predit = tracks(t).X_k_estime(1:3); 
        
        % --- ASSOCIATION MAHALANOBIS (CORRIGÉE POUR LE BPF) ---
        rayon_gating_max = [0.3; 0.3; 0.3]; % Rayon plafonné SPÉCIFIQUEMENT pour la porte d'association
        centers = zeros(3, N_boites);
        rayons = zeros(3, N_boites);
        for b = 1:N_boites
            centers(:, b) = mid(tracks(t).X_boxes(1:3, b));
            rayons(:, b) = min(rad(tracks(t).X_boxes(1:3, b)), rayon_gating_max); % On extrait la demi-largeur, plafonnée
        end
        
        % L'incertitude est la dispersion des centres + la largeur des boîtes (plafonnée)
        P_nuage = cov(centers') + diag(mean(rayons, 2).^2);
        
        S_k = P_nuage + R_cov;
        inv_S = inv(S_k);
        
        distances_mahal = zeros(nb_objets, 1);
        for d = 1:nb_objets
            if detections_associees(d)
                distances_mahal(d) = Inf;
            else
                innov_assoc = matrice_centres(d, :)' - centre_predit;
                distances_mahal(d) = sqrt(innov_assoc' * inv_S * innov_assoc);
            end
        end
        
        [min_dist, idx_best] = min(distances_mahal);
        
        if min_dist < 3.5
            mesure_Y = matrice_centres(idx_best, :)';
            Boite_LiDAR = midrad(mesure_Y, erreur_lidar_max);
            
            poids_temp = zeros(1, N_boites);
            
            for b = 1:N_boites
                Z_box = C_int * tracks(t).X_boxes(:, b);
                R_box = intersect(Z_box, Boite_LiDAR); 
                
                if any(isnan(R_box))
                    poids_temp(b) = 0;
                else
                    volume_Z = prod(rad(Z_box) * 2 + 1e-6);
                    volume_R = prod(rad(R_box) * 2 + 1e-6);
                    vraisemblance = volume_R / volume_Z;
                    
                    poids_temp(b) = tracks(t).poids(b) * vraisemblance;
                    
                    % 1. Contraction Position
                    tracks(t).X_boxes(1:3, b) = R_box;
                    
                    % 2. CONTRACTION VITESSE BACKWARD (LA RIGUEUR ENSEMBLISTE)
                    V_backward = (R_box - tracks(t).X_prec_estime(1:3)) / dt;
                    V_contractee = intersect(tracks(t).X_boxes(4:6, b), V_backward);
                    
                    if any(isnan(V_contractee))
                        tracks(t).X_boxes(4:6, b) = V_backward;
                    else
                        tracks(t).X_boxes(4:6, b) = V_contractee;
                    end
                end
            end
            
            somme_poids = sum(poids_temp);
            
            if somme_poids > 0
                tracks(t).poids = poids_temp / somme_poids;
                tracks(t).lost_frames = 0;
                
                N_eff = 1 / sum(tracks(t).poids.^2);
                if N_eff < N_boites / 2
                    tracks(t).X_boxes = subdivision_resampling(tracks(t).X_boxes, tracks(t).poids, N_boites);
                    tracks(t).poids = ones(1, N_boites) / N_boites;
                end
            else
                % SÉCURITÉ ANTI-FANTÔMES (la porte a déjà validé l'association, min_dist < 3.5)
                for b = 1:N_boites
                    tracks(t).X_boxes(1:3, b) = Boite_LiDAR;
                    V_backward = (Boite_LiDAR - tracks(t).X_prec_estime(1:3)) / dt;
                    tracks(t).X_boxes(4:6, b) = V_backward;
                end
                tracks(t).poids = ones(1, N_boites) / N_boites;
                tracks(t).lost_frames = 0;
            end
            
            centers = zeros(9, N_boites);
            for b = 1:N_boites
                centers(:, b) = mid(tracks(t).X_boxes(:, b));
            end
            tracks(t).X_k_estime = sum(centers .* repmat(tracks(t).poids, 9, 1), 2);
            
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k_estime(1:3)'];
            tracks(t).historique_observe = [tracks(t).historique_observe; mesure_Y'];
            tracks(t).age = tracks(t).age + 1;
            
            if tracks(t).age >= 10
                if ~tracks(t).confirmed 
                    tracks(t).display_id = next_display_id; 
                    next_display_id = next_display_id + 1;  
                end
                tracks(t).confirmed = true;
            end  
            detections_associees(idx_best) = true; 
            
        else
            tracks(t).lost_frames = tracks(t).lost_frames + 1;
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k_estime(1:3)'];
            tracks(t).historique_observe = [tracks(t).historique_observe; NaN, NaN, NaN]; 
        end
    end

    if ~isempty(tracks)
        pistes_valides = [tracks.lost_frames] <= 10;
        tracks = tracks(pistes_valides);
    end

    for d = 1:nb_objets
        if ~detections_associees(d)
            centre_initial = [matrice_centres(d, :)'; zeros(6,1)];
            incertitude = [erreur_lidar_max; erreur_lidar_max; erreur_lidar_max; 0.5; 0.5; 0.5; 0.1; 0.1; 0.1];
            boite_mere = midrad(centre_initial, incertitude);
            
            boites_initiales = diviser_boite(boite_mere, N_boites, 4); 
            poids_initiaux = ones(1, N_boites) / N_boites;
            
            nouvelle_piste = struct( ...
                'id', next_id, ...
                'X_boxes', boites_initiales, ...
                'X_prec_estime', centre_initial, ...
                'poids', poids_initiaux, ...
                'X_k_estime', centre_initial, ...
                'historique', matrice_centres(d, :), ...
                'historique_observe', matrice_centres(d, :), ...
                'lost_frames', 0, ...
                'age', 1, ...
                'confirmed', false, ...
                'display_id', 0 ...
            );
            tracks = [tracks, nouvelle_piste];
            next_id = next_id + 1; 
        end
    end
end

function new_boxes = subdivision_resampling(boxes, weights, N)
    cum_w = cumsum(weights);
    step = 1/N;
    u = rand() * step;
    counts = zeros(1, N);
    idx = 1;
    for i = 1:N
        while u > cum_w(idx) && idx < N
            idx = idx + 1;
        end
        counts(idx) = counts(idx) + 1;
        u = u + step;
    end
    
    new_boxes = intval(zeros(9, N));
    new_idx = 1;
    for i = 1:N
        k = counts(i);
        if k > 0
            rayons = rad(boxes(:, i));
            [~, dim_max] = max(rayons);
            
            sub_b = diviser_boite(boxes(:, i), k, dim_max);
            new_boxes(:, new_idx : new_idx+k-1) = sub_b;
            new_idx = new_idx + k;
        end
    end
end

function sub_boxes = diviser_boite(boite, k, dim)
    sub_boxes = repmat(boite, 1, k);
    if k == 1, return; end
    
    inf_val = inf(boite(dim));
    sup_val = sup(boite(dim));
    step = (sup_val - inf_val) / k;
    
    for i = 1:k
        new_inf = inf_val + (i-1)*step;
        new_sup = inf_val + i*step;
        sub_boxes(dim, i) = infsup(new_inf, new_sup);
    end
end