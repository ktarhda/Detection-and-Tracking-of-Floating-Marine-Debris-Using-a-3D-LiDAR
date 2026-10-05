function [tracks, detections_associees, next_display_id, next_id] = kalman_classique(tracks, matrice_centres, dt, A, C, Q, R, next_display_id, next_id)
    
    nb_objets = size(matrice_centres, 1);
    detections_associees = false(1, nb_objets);

    % 1. PRÉDICTION
    for t = 1:numel(tracks)
        tracks(t).X_k = A * tracks(t).X_k;
        tracks(t).P_k = A * tracks(t).P_k * A' + Q;
    end

    % 2. ASSOCIATION ET CORRECTION (Mahalanobis)
    for t = 1:numel(tracks)
        if nb_objets == 0, break; end
        
        S = C * tracks(t).P_k * C' + R;
        invS = inv(S); 
        
        diff = matrice_centres - tracks(t).X_k(1:3)';
        distances_cibles = sqrt(sum((diff * invS) .* diff, 2));
        
        distances_cibles(detections_associees) = Inf; 
        [min_dist, idx_best] = min(distances_cibles);
        
        if min_dist < 3.5 
            % CORRECTION
            mesure_Y = matrice_centres(idx_best, :)';
            Innovation = mesure_Y - (C * tracks(t).X_k);
            K = tracks(t).P_k * C' / S;
            
            tracks(t).X_k = tracks(t).X_k + K * Innovation;
            tracks(t).P_k = (eye(9) - K * C) * tracks(t).P_k;
            
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
            % OCCLUSION
            tracks(t).lost_frames = tracks(t).lost_frames + 1;
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k(1:3)'];
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
            nouvelle_piste = struct( ...
                'id', next_id, ...
                'X_k', [matrice_centres(d, :)'; 0; 0; 0; 0; 0; 0], ...
                'P_k', eye(9) * 10, ...
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



% Assosiation avec distance euclidienne 
%{
%Assosiation 
for t = 1:numel(tracks)
        if nb_objets == 0, break; end
        %{
   il prend la prédiction (tracks(t).X_k(1:3), là où il pense que l'objet est et soustrait les coordonnées de TOUTES les boîtes (matrice_centres).
 La ligne sqrt(sum(diff.^2, 2)) c'est le théorème de Pythagore en 3D pour calculer la distance euclidienne
        %}
        diff = matrice_centres - tracks(t).X_k(1:3)';
        distances_cibles = sqrt(sum(diff.^2, 2));
        distances_cibles(detections_associees) = Inf;
        [min_dist, idx_best] = min(distances_cibles);
        
        if min_dist < 2.5 %Si la distance entre la prédiction et la vraie boîte verte est inférieure à 2.5 mètres, l'algorithme accepte l'association. Si c'est plus grand, il refuse (l'objet ne peut pas se téléporter de 2.5m en 0.1 seconde)
            % CORRECTION DE KALMAN 
            mesure_Y = matrice_centres(idx_best, :)';
            Innovation = mesure_Y - (C * tracks(t).X_k);
            S = C * tracks(t).P_k * C' + R;
            K = tracks(t).P_k * C' / S;
            
            tracks(t).X_k = tracks(t).X_k + K * Innovation;
            tracks(t).P_k = (eye(6) - K * C) * tracks(t).P_k;
            
            % ENREGISTREMENT DES TRAJECTOIRES 
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k(1:3)'];
            % ON STOCKE LA MESURE BRUTE DU CAPTEUR (LIGNE JAUNE)
            tracks(t).historique_observe = [tracks(t).historique_observe; mesure_Y'];
            
            tracks(t).lost_frames = 0; % Piste confirmée visible
            
            % LOGIQUE DE CONFIRMATION AVEC LE DISPLAY_ID
            tracks(t).age = tracks(t).age + 1; % L'objet vieillit d'une frame
            if tracks(t).age >= 1
                if ~tracks(t).confirmed % Si c'est la toute première fois qu'il est confirmé
                    tracks(t).display_id = next_display_id; % On lui donne un beau numéro
                    next_display_id = next_display_id + 1;  % Le prochain aura le numéro suivant
                end
                tracks(t).confirmed = true;
            end  
            detections_associees(idx_best) = true; % Détection consommée
        else
            %Si aucune boîte n'était à moins de 2.5m, la cible est masquée. L'algorithme augmente le compteur lost_frames
            % L'objet est masqué : Kalman prédit, mais l'observation s'arrête (coupure de la ligne jaune)
            tracks(t).lost_frames = tracks(t).lost_frames + 1;
            tracks(t).historique = [tracks(t).historique; tracks(t).X_k(1:3)'];
            tracks(t).historique_observe = [tracks(t).historique_observe; NaN, NaN, NaN]; 
        end
    end

%}



