function [tracks, detections_associees, next_display_id, next_id] = filtre_particule(tracks, matrice_centres, dt, N_particules, Q_cov, R_cov, next_display_id, next_id)
    
    nb_objets = size(matrice_centres, 1);
    detections_associees = false(1, nb_objets);

   % Matrices du modèle cinématique à 9 ÉTATS
    A = [eye(3), eye(3)*dt, eye(3)*(0.5*dt^2);
         zeros(3,3), eye(3), eye(3)*dt;
         zeros(3,3), zeros(3,3), eye(3)];
    C = [eye(3), zeros(3,6)];
    
    % Inversion de R pour la correction des poids
    inv_R = inv(R_cov);

       % 1. PRÉDICTION (Déplacement de toutes les particules)
    ecart_max = [3; 3; 3; 2; 2; 2; 2; 2; 2]; % Dispersion max autorisée (position, vitesse, accélération)
    for t = 1:numel(tracks)
        % Ajout d'un bruit gaussien basé sur Q_cov
        bruit = (chol(Q_cov)' * randn(9, N_particules)); 
        
        % Déplacement selon le modèle physique + le bruit
        tracks(t).particules = A * tracks(t).particules + bruit;
        
        % L'état actuel estimé est la moyenne pondérée de toutes les particules
        tracks(t).X_k = sum(tracks(t).particules .* tracks(t).poids, 2);
        
        % SATURATION DE LA DISPERSION (évite l'explosion du nuage lors d'une occlusion)
        ecarts = tracks(t).particules - tracks(t).X_k;
        ecarts = max(min(ecarts, ecart_max), -ecart_max);
        tracks(t).particules = tracks(t).X_k + ecarts;
    end

    % 2. ASSOCIATION (MAHALANOBIS) ET CORRECTION
    for t = 1:numel(tracks)
        if nb_objets == 0, break; end
        
        centre_predit = tracks(t).X_k(1:3); % 3x1
        
        % --- CALCUL DE LA COVARIANCE DU NUAGE POUR MAHALANOBIS ---
        % On regarde comment les particules sont étalées (Incertitude de prédiction P)
        P_nuage = cov(tracks(t).particules(1:3, :)'); 
        
        % L'incertitude totale (S) = Incertitude du nuage + Incertitude du capteur
        S_covariance = P_nuage + R_cov; 
        inv_S = inv(S_covariance);
        
        % Calcul de la distance de Mahalanobis pour toutes les détections
        distances_mahal = zeros(nb_objets, 1);
        for d = 1:nb_objets
            if detections_associees(d)
                distances_mahal(d) = Inf; % Déjà associé à une autre piste
            else
                innov_assoc = matrice_centres(d, :)' - centre_predit;
                % Formule de la Distance de Mahalanobis
                distances_mahal(d) = sqrt(innov_assoc' * inv_S * innov_assoc);
            end
        end
        
        [min_mahal, idx_best] = min(distances_mahal);
        
        % --- SEUIL DE MAHALANOBIS ---
        % En statistiques (loi du Chi-2 à 3 degrés de liberté), une distance 
        % de Mahalanobis < 4.0 englobe plus de 99% des possibilités physiques réelles.
        if min_mahal < 3.5
            mesure_Y = matrice_centres(idx_best, :)';
            
            % --- CORRECTION DES POIDS ---
            % On calcule l'écart entre chaque particule et la mesure LiDAR
            innov_particules = mesure_Y - C * tracks(t).particules; % 3 x N_particules
            
            % Calcul de la distance de Mahalanobis pour chaque particule
            dist2_particules = sum(innov_particules .* (inv_R * innov_particules), 1); 
            
            % Calcul de la Vraisemblance (Likelihood) gaussienne
            vraisemblance = exp(-0.5 * dist2_particules);
            vraisemblance = vraisemblance + 1e-15; % Sécurité division par zéro
            
            % Mise à jour des poids
            tracks(t).poids = tracks(t).poids .* vraisemblance;
            tracks(t).poids = tracks(t).poids / sum(tracks(t).poids); % Normalisation
            
            % --- RÉÉCHANTILLONNAGE (Resampling) ---
            N_eff = 1 / sum(tracks(t).poids.^2);
            if N_eff < N_particules / 2
                indices = randsample(1:N_particules, N_particules, true, tracks(t).poids);
                tracks(t).particules = tracks(t).particules(:, indices);
                tracks(t).poids = ones(1, N_particules) / N_particules; 
            end
            
            % Estimation finale après correction
            tracks(t).X_k = sum(tracks(t).particules .* tracks(t).poids, 2);
            
            % Sauvegardes
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k(1:3)'];
            tracks(t).historique_observe = [tracks(t).historique_observe; mesure_Y'];
            tracks(t).lost_frames = 0; 
            tracks(t).age = tracks(t).age + 1;
            
            if tracks(t).age >= 5
                if ~tracks(t).confirmed 
                    tracks(t).display_id = next_display_id; 
                    next_display_id = next_display_id + 1;  
                end
                tracks(t).confirmed = true;
            end  
            detections_associees(idx_best) = true; 
            
        else
            % --- COASTING (OBJET PERDU) ---
            tracks(t).lost_frames = tracks(t).lost_frames + 1;
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k(1:3)'];
            tracks(t).historique_observe = [tracks(t).historique_observe; NaN, NaN, NaN]; 
        end
    end

    % 3. NETTOYAGE DES FANTÔMES
    if ~isempty(tracks)
        pistes_valides = [tracks.lost_frames] <= 10;
        tracks = tracks(pistes_valides);
    end

    % 4. CRÉATION D'UNE NOUVELLE PISTE
    for d = 1:nb_objets
        if ~detections_associees(d)
            centre_initial = [matrice_centres(d, :)'; zeros(6,1)];
            % Matrice d'incertitude initiale
            P_init = diag([0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1].^2); 
            nouvelles_particules = centre_initial + (chol(P_init)' * randn(9, N_particules));

            nouveaux_poids = ones(1, N_particules) / N_particules;
            
            nouvelle_piste = struct( ...
                'id', next_id, ...
                'X_k', centre_initial, ...
                'particules', nouvelles_particules, ...
                'poids', nouveaux_poids, ...
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