clear; clc; close all;
disp('  PERCEPTION MARINE : DÉTECTION  ');
 
%pcapFile='un_seau.pcap'; jsonFile='un_seau.json';
%pcapFile='un_seau_version2.pcap'; jsonFile='un_seau_version2.json';
pcapFile='le_boit.pcap'; jsonFile='le_boit.json';
%pcapFile='3objets.pcap'; jsonFile='3objets.json';
% LECTURE JSON + RÉSOLUTION
config = jsondecode(fileread(jsonFile));

if isfield(config, 'data_format')
    N_cols   = config.data_format.columns_per_frame;
    N_rows   = config.data_format.pixels_per_column;
    phi_brut = config.beam_altitude_angles;
elseif isfield(config, 'lidar_data_format')
    N_cols   = config.lidar_data_format.columns_per_frame;
    N_rows   = config.lidar_data_format.pixels_per_column;
    phi_brut = config.beam_intrinsics.beam_altitude_angles;
else
    error('Format JSON non reconnu. Vérifiez le fichier de configuration du LiDAR.');
end

theta_brut = (0 : N_cols-1) * (360 / N_cols);
% ROI 
demi_angle  = 90; % on peux l'augmenter de  5° à 90°
theta_avant = 180;
ROI         = abs(theta_brut - theta_avant) <= demi_angle;
col_ROI     = find(ROI);

% LECTEUR + IMU
reader   = ousterFileReader(pcapFile, jsonFile);
imu_data = readIMU(reader);
total_frames = reader.NumberOfFrames;

% VISUALISATEUR 
viewer = pcplayer([-30 30], [-30 30], [-5 5]);
title(viewer.Axes, 'Détection Objets Aquatiques');

% FILTRE DE KALMAN MULTI-OBJETS (MOT) - INITIALISATION
mode_str = ''; 
if isfield(config, 'config_params') && isfield(config.config_params, 'lidar_mode')
    mode_str = config.config_params.lidar_mode; 
elseif isfield(config, 'lidar_mode')
    mode_str = config.lidar_mode; 
end
if ~isempty(mode_str)
    parties = split(mode_str, 'x'); 
    frequence_hz = str2double(parties{2}); 
    dt = 1.0 / frequence_hz; 
else
    % Sécurité si le fichier est corrompu ou illisible
    disp('Attention: lidar_mode introuvable. Fréquence par défaut 10 Hz.');
    frequence_hz = 10;
    dt = 0.1;
end

fprintf('-> Fréquence LiDAR détectée : %d Hz | dt = %.3f seconde\n', frequence_hz, dt);

A = [eye(3), eye(3)*dt; 
     zeros(3,3), eye(3)];

C = [eye(3), zeros(3,3)];

Q = eye(6) * 0.01; 

R = eye(3) * 0.03; % Bruit du capteur

% Structure dynamique enrichie avec la logique de Confirmation (Age)
tracks = struct('id', {}, 'X_k', {}, 'P_k', {}, 'historique', {}, ...
                'historique_observe', {}, 'lost_frames', {}, ...
                'age', {}, 'confirmed', {});
next_id = 1; % Compteur pour donner des numéros (ID 1, ID 2, ID 3...)

% BOUCLE PRINCIPALE
frame_idx = 0;
for i = 1 : (total_frames -334)
    if ~isOpen(viewer)
        break;
    end
    frame_idx = frame_idx + 1;
    % 1. ACQUISITION + CALIBRATION
    [ptCloud_brut, pcAttributes] = readFrame(reader, i);
    rho_brut       = double(pcAttributes.Range);
    ptCloud_aligne = correct_inclinaison(ptCloud_brut, imu_data);
    % 2. APPLICATION DU ROI 
    X_aligne = ptCloud_aligne.Location(:, col_ROI, 1);
    Y_aligne = ptCloud_aligne.Location(:, col_ROI, 2);
    Z_aligne = ptCloud_aligne.Location(:, col_ROI, 3); 
    % Le Rho de travail est directement le Rho brut de la zone d'intérêt
    portee_min = -0.4;
    portee_max = -15.0;
    limite_gauche = -7.0;
    rho_sans_sol = rho_brut(:, col_ROI); 
    zone_lac = (X_aligne >= portee_max & X_aligne <= portee_min) & (Y_aligne >= limite_gauche);
    rho_sans_sol(~zone_lac) = 0;
   % Le nuage d'étude est directement la matrice du ROI
    matrice_3D_ROI = cat(3, X_aligne, Y_aligne, Z_aligne);
    ptCloud_sans_sol = pointCloud(matrice_3D_ROI);
   
    % 3. SEGMENTATION DIRECTE SUR RHO BRUT
    [candidats, boites] = segmentation_rho_marin(rho_sans_sol, ptCloud_sans_sol);
    nb_objets = numel(candidats);
    fprintf('Objets détectés sur la frame %d : %d\n', frame_idx, nb_objets);
    
    % 4. COLORATION (FOND SELON Rang + OBJETS EN BLEU)
    ptCloud_visu = ptCloud_aligne; 
    [M, N, ~] = size(ptCloud_visu.Location);
    %  4.1. CRÉATION DU FOND COLORÉ (SELON LA DISTANCE RHO) 
    % On limite l'échelle de couleur à 40 mètres pour avoir un bon contraste
    distance_max_couleur = 40.0; 
    rho_norm = rho_brut / distance_max_couleur;
    rho_norm(rho_norm > 1) = 1; % On bloque le maximum à 1 (tout ce qui est > 40m aura la même couleur max)
    % On choisit une palette de couleurs MATLAB ('jet', 'parula', 'turbo' ou 'hsv')
    % 'jet' donne un bel arc-en-ciel du bleu (proche) au rouge (loin).
    palette = jet(256); 
    % On convertit nos distances en indices pour la palette (de 1 à 256)
    indices_couleurs = round(rho_norm * 255) + 1;
    % On construit l'image RGB de fond
    colors = zeros(M, N, 3, 'uint8');
    colors(:,:,1) = reshape(palette(indices_couleurs, 1) * 255, M, N);
    colors(:,:,2) = reshape(palette(indices_couleurs, 2) * 255, M, N);
    colors(:,:,3) = reshape(palette(indices_couleurs, 3) * 255, M, N);
    %  4.2. SUPERPOSITION DES OBJETS EN BLEU VIF
    X_visu = ptCloud_visu.Location(:,:,1);
    Y_visu = ptCloud_visu.Location(:,:,2);
    Z_visu = ptCloud_visu.Location(:,:,3);
    
    if ~isempty(boites)
        for b = 1:size(boites, 1)
            box = boites(b, :);
            xMin = box(1) - box(4)/2; xMax = box(1) + box(4)/2;
            yMin = box(2) - box(5)/2; yMax = box(2) + box(5)/2;
            zMin = box(3) - box(6)/2; zMax = box(3) + box(6)/2;
            
            % Masque : 1 si le point appartient à la boîte
            masque_boite = (X_visu >= xMin & X_visu <= xMax) & ...
                           (Y_visu >= yMin & Y_visu <= yMax) & ...
                           (Z_visu >= zMin & Z_visu <= zMax);
            
            % On écrase le fond arc-en-ciel pour mettre les objets en bleu pur
            canal_R = colors(:,:,1); canal_G = colors(:,:,2); canal_B = colors(:,:,3);
            
            canal_R(masque_boite) = 0;   % R = 0
            canal_G(masque_boite) = 0;   % G = 0
            canal_B(masque_boite) = 255; % B = 255
            
            colors(:,:,1) = canal_R; colors(:,:,2) = canal_G; colors(:,:,3) = canal_B;
        end
    end
    % On applique la couleur finale
    ptCloud_visu.Color = colors;
    % 4.5. ALGORITHME DE MULTI-OBJECT TRACKING (MOT)
    % 1. PRÉDICTION
    for t = 1:numel(tracks)
        tracks(t).X_k = A * tracks(t).X_k;
        tracks(t).P_k = A * tracks(t).P_k * A' + Q;
    end
%{
on utilise reshape pour extraire uniquement les centres (X, Y, Z) et en faire un tableau mathématique simple (une matrice).
detections_associees : C'est un tableau rempli de "Faux" (false).
 Il sert de "feuille d'émargement". À chaque fois qu'une boîte verte sera attribuée à une trajectoire, on cochera "Vrai". 
Cela empêche que deux trajectoires différentes s'accrochent à la même boîte verte.
%}
    matrice_centres = [];
    if nb_objets > 0
        matrice_centres = reshape([candidats.centre], 3, [])';
    end
    detections_associees = false(1, nb_objets);

    % 2. ASSOCIATION
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
            
            % LOGIQUE DE CONFIRMATION 
            tracks(t).age = tracks(t).age + 1; % L'objet vieillit d'une frame
            if tracks(t).age >= 10
                % Si l'objet a survécu 10 frames consécutives, c'est un VRAI objet !
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
    % 3. NETTOYAGE
    if ~isempty(tracks)
        pistes_valides = [tracks.lost_frames] <= 10;
        tracks = tracks(pistes_valides);
    end
   % 4. CRÉATION D'UNE NOUVELLE PISTE
    for d = 1:nb_objets
        if ~detections_associees(d)
            nouvelle_piste.id = next_id;
            nouvelle_piste.X_k = [matrice_centres(d, :)'; 0; 0; 0]; %0.0.0 au debut l'Algo ne connais pas les vitesse
            nouvelle_piste.P_k = eye(6) * 10;
            nouvelle_piste.historique = matrice_centres(d, :);
            nouvelle_piste.historique_observe = matrice_centres(d, :);
            nouvelle_piste.lost_frames = 0;
            
            % INITIALISATION DE L'ÂGE 
            nouvelle_piste.age = 1;
            nouvelle_piste.confirmed = false; % Piste "Tentative" (Cachée)
            tracks = [tracks, nouvelle_piste];
            next_id = next_id + 1; 
        end
    end
    % 5. RECADRAGE STRICT 
    X_final = ptCloud_visu.Location(:,:,1);
    Y_final = ptCloud_visu.Location(:,:,2);
    Z_final = ptCloud_visu.Location(:,:,3);
    % On recadre proprement à 50 mètres pour matcher avec le pcplayer
    idx_valides = find( X_final >= -50 & X_final <= 50 & ...
                        Y_final >= -50 & Y_final <= 50 & ...
                        Z_final >= -5   & Z_final <= 5 );
                        
    ptCloud_visu = select(ptCloud_visu, idx_valides); 
    % 6. AFFICHAGE 3D + DESSIN DES BOÎTES VERTES
    view(viewer, ptCloud_visu);
    delete(findobj(viewer.Axes, 'Type', 'Patch')); % Nettoyage anciennes boîtes

    if ~isempty(boites)
        couleurs = repmat([0 1 0], size(boites, 1), 1); % Vert éclatant
        showShape('cuboid', boites, ...
                  'Parent', viewer.Axes, ...
                  'Color', couleurs, ...
                  'Opacity', 0.1, ...
                  'LineWidth', 0.3);
    end
    % AFFICHAGE DES TRAJECTOIRES MULTIPLES (MOT)
    % Nettoyage complet
    delete(findobj(viewer.Axes, 'Tag', 'LigneTrajectoireMOT'));
    delete(findobj(viewer.Axes, 'Tag', 'LigneObserveeMOT'));
    delete(findobj(viewer.Axes, 'Tag', 'PointActuelMOT'));
    delete(findobj(viewer.Axes, 'Tag', 'TexteIDMOT'));

    if ~isempty(tracks)
        hold(viewer.Axes, 'on');
          
       for t = 1:numel(tracks)
            %  ON NE DESSINE QUE LES VRAIS OBJETS 
            if tracks(t).confirmed
                
                hist_filtre  = tracks(t).historique;
                hist_observe = tracks(t).historique_observe;
                
                % Calcul de la vitesse
                vitesse_m_s = norm(tracks(t).X_k(4:6));
                
                % Dessin de la ligne rouge
                if size(hist_filtre, 1) > 1
                    plot3(viewer.Axes, hist_filtre(:,1), hist_filtre(:,2), hist_filtre(:,3), ...
                          'r-', 'LineWidth', 2.5, 'Tag', 'LigneTrajectoireMOT');
                end
                
                % Dessin de la ligne jaune observée
                if size(hist_observe, 1) > 1
                    plot3(viewer.Axes, hist_observe(:,1), hist_observe(:,2), hist_observe(:,3), ...
                          'y--', 'LineWidth', 1.5, 'Tag', 'LigneObserveeMOT');
                end
                
                % Point rouge actuel
                plot3(viewer.Axes, tracks(t).X_k(1), tracks(t).X_k(2), tracks(t).X_k(3), ...
                      'ro', 'MarkerSize', 3.0, 'MarkerFaceColor', 'r', 'Tag', 'PointActuelMOT');
                      
                % Texte ID et vitesse
                texte_affichage = sprintf('ID %d | %.2f m/s', tracks(t).id, vitesse_m_s);
                text(viewer.Axes, tracks(t).X_k(1), tracks(t).X_k(2), tracks(t).X_k(3) + 0.6, ...
                     texte_affichage, 'Color', 'green', 'FontSize', 8, 'FontWeight', 'bold', 'Tag', 'TexteIDMOT');
                     
            end % <-- Fin du if tracks(t).confirmed
        end
        hold(viewer.Axes, 'off');
    end
    drawnow limitrate;
end
disp('Fin de l''acquisition marine.');