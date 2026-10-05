function [tracks, detections_associees, next_display_id, next_id] = kalman_ensembliste(tracks, matrice_centres, dt, A, C, Q, R, V_modele, erreur_lidar_max, next_display_id, next_id)
    
    nb_objets = size(matrice_centres, 1);
    detections_associees = false(1, nb_objets);

    % OPTIMISATION 1 : MATRICE Q CINÉMATIQUE OPTIMALE
    q_bruit = 0.01; % Densité du bruit (vagues)
    Q_opt = [ (dt^5)/20*eye(3), (dt^4)/8*eye(3),  (dt^3)/6*eye(3) ; 
              (dt^4)/8*eye(3),  (dt^3)/3*eye(3),  (dt^2)/2*eye(3) ; 
              (dt^3)/6*eye(3),  (dt^2)/2*eye(3),  dt*eye(3)       ] * q_bruit;

    % Les matrices du modèle d'état sont en intervalles
    A_int = intval(A);
    C_int = intval(C);
    
    % (Note : Q_opt et R restent des matrices réelles pour le calcul de la borne P+
    % afin d'éviter l'explosion des intervalles par "Wrapping Effect". 
    % Leurs équivalents en intervalles sont déjà V_modele et erreur_lidar_max)

    % Initialisation de la borne P+ et de l'historique
    if ~isfield(tracks, 'P_k_plus')
        if isempty(tracks)
            tracks = struct('id', {}, 'X_k_intval', {}, 'P_k_plus', {}, 'X_prec_intval', {}, ...
                            'historique', {}, 'historique_observe', {}, 'lost_frames', {}, ...
                            'age', {}, 'confirmed', {}, 'display_id', {});
        else
            for j = 1:numel(tracks)
                tracks(j).P_k_plus = eye(9); 
                tracks(j).X_prec_intval = tracks(j).X_k_intval;
            end
        end
    end

    % 1. PRÉDICTION ENSEMBLISTE
    for t = 1:numel(tracks)
        tracks(t).X_prec_intval = tracks(t).X_k_intval; % Sauvegarde pour le Forward-Backward
        
        % Prédiction de l'état strict en intervalles
        tracks(t).X_k_intval = A_int * tracks(t).X_k_intval + V_modele;
        
        % Prédiction de la borne d'incertitude supérieure (Matrice réelle)
        tracks(t).P_k_plus = A * tracks(t).P_k_plus * A' + Q_opt; 
    end

    % 2. ASSOCIATION (MAHALANOBIS) ET CORRECTION
    for t = 1:numel(tracks)
        if nb_objets == 0, break; end
        
        % Centre de la boîte prédite
        centre_predit = mid(tracks(t).X_k_intval(1:3)); 
        
        % --- CALCUL DE LA MATRICE S POUR MAHALANOBIS ---
        % On utilise l'incertitude théorique P+ pour l'association
        S_k = C * tracks(t).P_k_plus * C' + R;
        inv_S = inv(S_k);
        
        distances_mahal = zeros(nb_objets, 1);
        
        for d = 1:nb_objets
            if detections_associees(d)
                distances_mahal(d) = Inf; % Déjà associé
            else
                innov_assoc = matrice_centres(d, :)' - centre_predit;
                % Distance de Mahalanobis
                distances_mahal(d) = sqrt(innov_assoc' * inv_S * innov_assoc);
            end
        end
        
        [min_mahal, idx_best] = min(distances_mahal);
        
        % Seuil statistique de Mahalanobis à 4.0 (Loi du Chi-2 à 3 ddl = >99% de confiance)
        if min_mahal < 3.5
            mesure_Y = matrice_centres(idx_best, :)';
            
            % Création de la mesure en boîte d'intervalle
            Boite_LiDAR = midrad(mesure_Y, erreur_lidar_max);

            % L'ALGORITHME UBIKF & CONTRACTEUR (Xiong 2013 / Tran 2021)
            
            % 1. Calcul du Gain K optimal à partir de la borne supérieure P+
            K = tracks(t).P_k_plus * C' / S_k;
            
            % 2. Les équations a1, a2, a3 EXACTES de l'article (Xiong et al. 2013)
            a1 = C_int * tracks(t).X_k_intval;
            a2 = Boite_LiDAR - a1;
            a3 = K * a2;
            X_temp = tracks(t).X_k_intval + a3;
            
            % 3. Mise à jour de la borne supérieure P+ (Équation 22 de Tran et al. 2021)
            tracks(t).P_k_plus = (eye(9) - K * C) * tracks(t).P_k_plus;

            % OPTIMISATION 2 : PROPAGATION DE CONTRAINTES (Xiong et al. 2013)
            Nouvelle_Position = intersect(X_temp(1:3), Boite_LiDAR);
            if any(isnan(Nouvelle_Position)) 
                Nouvelle_Position = Boite_LiDAR;
            end
            tracks(t).X_k_intval(1:3) = Nouvelle_Position;
            
            % Contraction de la vitesse (Backward Kinematic Constraint)
            V_backward = (Nouvelle_Position - tracks(t).X_prec_intval(1:3)) / dt;
         
            % OPTIMISATION 3 : LARGEUR MATHÉMATIQUE (1.5 Sigma = 86.6%)
            vitesse_mise_a_jour = mid(X_temp(4:6));
            sigma_vx = sqrt(tracks(t).P_k_plus(4,4));
            sigma_vy = sqrt(tracks(t).P_k_plus(5,5));
            sigma_vz = sqrt(tracks(t).P_k_plus(6,6));
            
            incertitude_V = 1.5 * [sigma_vx; sigma_vy; sigma_vz]; 
            V_ubikf = midrad(vitesse_mise_a_jour, incertitude_V);
            
            % L'intersection absolue : Théorie UBIKF combinée au Forward-Backward
            Vitesse_Contractee = intersect(V_ubikf, V_backward);
            
            if any(isnan(Vitesse_Contractee))
                tracks(t).X_k_intval(4:6) = V_ubikf;
            else
                tracks(t).X_k_intval(4:6) = Vitesse_Contractee;
            end
            % --- NOUVEAU : CONTRAINTE ACCÉLÉRATION ---
            % 1. Calcul de l'accélération backward
            A_backward = (tracks(t).X_k_intval(4:6) - tracks(t).X_prec_intval(4:6)) / dt;
            
            % 2. Extraction du centre mathématique de Kalman pour l'accélération
            accel_mise_a_jour = mid(X_temp(7:9));
            
            % 3. Extraction de l'incertitude théorique garantie par P+
            sigma_ax = sqrt(tracks(t).P_k_plus(7,7));
            sigma_ay = sqrt(tracks(t).P_k_plus(8,8));
            sigma_az = sqrt(tracks(t).P_k_plus(9,9));
            incertitude_A = 1.5 * [sigma_ax; sigma_ay; sigma_az]; 
            
            A_ubikf = midrad(accel_mise_a_jour, incertitude_A);
            
            % 4. Intersection absolue
            Accel_Contractee = intersect(A_ubikf, A_backward);
            if any(isnan(Accel_Contractee))
                tracks(t).X_k_intval(7:9) = A_ubikf;
            else
                tracks(t).X_k_intval(7:9) = Accel_Contractee;
            end
            
            %  SAUVEGARDES
            tracks(t).historique = [tracks(t).historique; mid(tracks(t).X_k_intval(1:3))'];
            tracks(t).historique_observe = [tracks(t).historique_observe; mesure_Y'];
            
            tracks(t).lost_frames = 0; 
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
            % COASTING
            sigma_vx = sqrt(tracks(t).P_k_plus(4,4));
            sigma_vy = sqrt(tracks(t).P_k_plus(5,5));
            sigma_vz = sqrt(tracks(t).P_k_plus(6,6));
            
            incertitude_V = 1.5 * [sigma_vx; sigma_vy; sigma_vz];
            vitesse_gelee = mid(tracks(t).X_k_intval(4:6));
            
            tracks(t).X_k_intval(4:6) = midrad(vitesse_gelee, incertitude_V);
            % --- NOUVEAU : COASTING DE L'ACCÉLÉRATION ---
            sigma_ax = sqrt(tracks(t).P_k_plus(7,7));
            sigma_ay = sqrt(tracks(t).P_k_plus(8,8));
            sigma_az = sqrt(tracks(t).P_k_plus(9,9));
            incertitude_A = 1.5 * [sigma_ax; sigma_ay; sigma_az];
            
            accel_gelee = mid(tracks(t).X_k_intval(7:9));
            tracks(t).X_k_intval(7:9) = midrad(accel_gelee, incertitude_A);
            tracks(t).lost_frames = tracks(t).lost_frames + 1;
            tracks(t).historique = [tracks(t).historique; mid(tracks(t).X_k_intval(1:3))'];
            tracks(t).historique_observe = [tracks(t).historique_observe; NaN, NaN, NaN]; 
        end
    end

    % 3. NETTOYAGE
    if ~isempty(tracks)
        pistes_valides = [tracks.lost_frames] <= 10;
        tracks = tracks(pistes_valides);
    end

    % 4. CRÉATION D'UNE NOUVELLE PISTE
    for d = 1:nb_objets
        if ~detections_associees(d)
            centre_initial = [matrice_centres(d, :)'; 0; 0; 0; 0; 0; 0];
            incertitude_initiale = [erreur_lidar_max; erreur_lidar_max; erreur_lidar_max; 0.5; 0.5; 0.5; 0.1; 0.1; 0.1];
            
            nouvelle_piste = struct( ...
                'id', next_id, ...
                'X_k_intval', midrad(centre_initial, incertitude_initiale), ...
                'P_k_plus', eye(9), ...
                'X_prec_intval', midrad(centre_initial, incertitude_initiale), ...
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