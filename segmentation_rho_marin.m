function [candidats, boites] = segmentation_rho_marin(rho_sans_sol, ptCloud_sans_sol)
   
    [lignes, cols] = size(rho_sans_sol);
    
    % On prépare des matrices vides (pleines de 0) pour nos frontières
    frontiere_horiz = false(lignes, cols-1);
    frontiere_vert  = false(lignes-1, cols);
    
    seuil_initial = 0.05; 
    precision_ouster = 0.03;

    % 1 BALAYAGE HORIZONTAL (Gauche à Droite) 
    for i = 1:lignes
        seuil_actuel = seuil_initial; % On réinitialise le seuil au début de chaque ligne
        historique_sauts = [];  
        for j = 2:cols
            if rho_sans_sol(i,j)> 0 && rho_sans_sol(i,j-1)> 0
                % On calcule le saut entre le pixel et son voisin de gauche
                saut = abs(rho_sans_sol(i, j) - rho_sans_sol(i, j-1));
                
                if isempty(historique_sauts)
                    condition = saut > seuil_actuel;
                else
                    % 1.5 pour mieux tolérer le bruit naturel du LiDAR
                    condition = saut > (1.5 * seuil_actuel); 
                end
                
                if condition
                    % SAUT DÉTECTÉ ! On met un 1 sur la frontière
                    frontiere_horiz(i, j-1) = true;
                    
                    % On réinitialise le seuil
                    seuil_actuel = seuil_initial; 
                    historique_sauts = [];
                else
                    % PAS DE SAUT ! Surface continue → on enrichit la moyenne
                    historique_sauts = [historique_sauts, saut];
                    seuil_actuel = max(precision_ouster,mean(historique_sauts));
                end
            end
        end
    end
    
    % 2 BALAYAGE VERTICAL (Haut en Bas) 
    for j = 1:cols
        seuil_actuel = seuil_initial; % On réinitialise au début de chaque colonne
        historique_sauts = []; 

        for i = 2:lignes
           if rho_sans_sol(i,j)> 0 && rho_sans_sol(i-1,j)> 0
                % On calcule le saut entre le pixel et son voisin du dessus
                saut = abs(rho_sans_sol(i, j) - rho_sans_sol(i-1, j));
                
                if isempty(historique_sauts)
                    condition = saut > seuil_actuel;
                else
                    condition = saut > (1.5 * seuil_actuel); 
                end
                
                if condition  
                    frontiere_vert(i-1, j) = true;
                    seuil_actuel = seuil_initial;
                    historique_sauts = [];
                else
                   historique_sauts = [historique_sauts, saut];
                   seuil_actuel     = max(precision_ouster,mean(historique_sauts));
                end
           end
        end
    end
     
    % Image des objets valides
    image_objets = rho_sans_sol > 0;
    
    % Le coup de ciseaux : on coupe là où il y a des frontières
    image_objets(:, 1:end-1) = image_objets(:, 1:end-1) & ~frontiere_horiz;
    image_objets(:, 2:end)   = image_objets(:, 2:end)   & ~frontiere_horiz;
    image_objets(1:end-1, :) = image_objets(1:end-1, :) & ~frontiere_vert;
    image_objets(2:end, :)   = image_objets(2:end, :)   & ~frontiere_vert;
 
   
    % LA "COLLE" MARITIME (Fermeture Morphologique)
    % Rebouche les petits trous causés par l'absorption de l'eau sur le plastique
    pinceau = strel('rectangle', [3 3]); 
    image_objets = imclose(image_objets, pinceau);
   
    
    % 2. ÉTIQUETAGE
    % On passe en 8 connexions (au lieu de 4) pour accrocher les pixels lointains en diagonale
    [labels_2D, nb_regions] = bwlabel(image_objets, 8);
    
    % 3. CRÉATION DES BOÎTES 
    X = ptCloud_sans_sol.Location(:, :, 1);
    Y = ptCloud_sans_sol.Location(:, :, 2);
    Z = ptCloud_sans_sol.Location(:, :, 3);
    
    candidats = [];
    boites    = [];
    
    for k = 1 : nb_regions
        masque_k = (labels_2D == k);
        X_k = X(masque_k); Y_k = Y(masque_k); Z_k = Z(masque_k);
        
        % Nettoyage des NaN
        ok = ~isnan(X_k) & ~isnan(Y_k) & ~isnan(Z_k);
        X_k = X_k(ok); Y_k = Y_k(ok); Z_k = Z_k(ok);
        
        % Filtre anti-bruit : On abaisse à 3 points min car les objets marins au loin renvoient peu de laser
        if length(X_k) < 5
            continue;
        end
        
        % On calcule les dimensions cartésiennes de l'objet
        larg = max(X_k) - min(X_k);
        prof = max(Y_k) - min(Y_k);
        haut = max(Z_k) - min(Z_k);
        
        % Le centre parfait (géométrique)
        x_centre = min(X_k) + (larg / 2);
        y_centre = min(Y_k) + (prof / 2);
        z_centre = min(Z_k) + (haut / 2);
        
        % Enregistrement pour le Kalman (avec le nb de points si besoin)
        candidat.centre = [x_centre, y_centre, z_centre];
        candidat.nb_imp = length(X_k);
        candidats = [candidats, candidat];
        
        % Création de la boîte verte pour l'affichage
        boites = [boites; x_centre, y_centre, z_centre, larg, prof, haut, 0, 0, 0];
    end

end