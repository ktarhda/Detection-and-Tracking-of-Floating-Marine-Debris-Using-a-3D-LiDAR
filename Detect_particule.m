clear; clc; close all;
disp('  PERCEPTION MARINE : FILTRES PARTICULAIRES (Classique & Boîtes) ');
fenetre_lissage = 9;   % Taille de la fenêtre du filtre moyenne glissante (en frames)
% =========================================================================
% SÉLECTION DU FILTRE
% =========================================================================
% Choisissez : 'PARTICULAIRE' (Nuage de points) ou 'BOX_PARTICULAIRE' (Boîtes INTLAB)
%CHOIX_FILTRE = 'BOX_PARTICULAIRE'; 
CHOIX_FILTRE = 'PARTICULAIRE';

if strcmp(CHOIX_FILTRE, 'BOX_PARTICULAIRE')
    disp('-> Démarrage INTLAB pour le Box Particle Filter...');
    startintlab; 
end

% PARAMÈTRES COMMUNS
dt = 0.1;
pcapFile='3objets.pcap'; jsonFile='3objets.json';
%pcapFile='un_seau.pcap'; jsonFile='un_seau.json';

% LECTURE JSON + RÉSOLUTION
config = jsondecode(fileread(jsonFile));
if isfield(config, 'data_format')
    N_cols   = config.data_format.columns_per_frame;
    N_rows   = config.data_format.pixels_per_column;
else
    N_cols   = config.lidar_data_format.columns_per_frame;
    N_rows   = config.lidar_data_format.pixels_per_column;
end
theta_brut = (0 : N_cols-1) * (360 / N_cols);

% ROI 
demi_angle  = 90; 
theta_avant = 180;
ROI         = abs(theta_brut - theta_avant) <= demi_angle;
col_ROI     = find(ROI);

% LECTEUR + IMU
reader   = ousterFileReader(pcapFile, jsonFile);
imu_data = readIMU(reader);
total_frames = reader.NumberOfFrames;

% VISUALISATEUR 
viewer = pcplayer([-30 30], [-30 30], [-5 5]);
%title(viewer.Axes, ['Détection Objets Aquatiques : ' CHOIX_FILTRE]);

% =========================================================================
% INITIALISATION SPÉCIFIQUE (9 ÉTATS)
% =========================================================================
if strcmp(CHOIX_FILTRE, 'PARTICULAIRE')
    N_particules = 1000; 
    % Bruit sur 9 dimensions (pos, vit, acc)
    Q_cov = diag([0.02, 0.02, 0.02, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05].^2); 
    R_cov = diag([0.03, 0.03, 0.03].^2);
    tracks = struct('id', {}, 'X_k', {}, 'particules', {}, 'poids', {}, 'historique', {}, ...
                    'historique_observe', {}, 'lost_frames', {}, ...
                    'age', {}, 'confirmed', {}, 'display_id', {});
elseif strcmp(CHOIX_FILTRE, 'BOX_PARTICULAIRE')
    N_boites = 10; 
    erreur_lidar_max = 0.03; 
    R_cov = diag([0.03, 0.03, 0.03].^2); 
    % Modèle INTLAB sur 9 dimensions
    V_modele = midrad(zeros(9,1), [0.02; 0.02; 0.02; 0.1; 0.1; 0.1; 0.5; 0.5; 0.5]);
    tracks = struct('id', {}, 'X_boxes', {}, 'X_prec_estime', {}, 'poids', {}, 'X_k_estime', {}, 'historique', {}, ...
                    'historique_observe', {}, 'lost_frames', {}, ...
                    'age', {}, 'confirmed', {}, 'display_id', {});
end

next_id = 1; 
next_display_id = 1; 

archive_vitesses = struct('display_id', {}, 'vitesses_moy', {}, 'vitesses_min', {}, 'vitesses_max', {}, ...
                          'vx_moy', {}, 'vx_min', {}, 'vx_max', {}, ...
                          'vy_moy', {}, 'vy_min', {}, 'vy_max', {}, ...
                          'vz_moy', {}, 'vz_min', {}, 'vz_max', {}, ...
                          'ax_moy', {}, 'ax_min', {}, 'ax_max', {}, ...
                          'ay_moy', {}, 'ay_min', {}, 'ay_max', {}, ...
                          'az_moy', {}, 'az_min', {}, 'az_max', {}, ...
                          'frames', {}, 'est_observe', {});

% BOUCLE PRINCIPALE
frame_idx = 0;
% MÉTRIQUE DE CARDINALITÉ
n_reel = 3;                  % Nombre réel de cibles déployées dans le canal
historique_cardinalite = []; % Nombre de pistes confirmées à chaque trame
for i = 720 : (total_frames - 470)
    if ~isOpen(viewer)
        break;
    end
    frame_idx = frame_idx + 1;
    
    % 1. ACQUISITION + CALIBRATION
    [ptCloud_brut, pcAttributes] = readFrame(reader, i);
    rho_brut = double(pcAttributes.Range);
    ptCloud_aligne = correct_inclinaison(ptCloud_brut, imu_data);
    
    % 2. APPLICATION DU ROI 
    X_aligne = ptCloud_aligne.Location(:, col_ROI, 1);
    Y_aligne = ptCloud_aligne.Location(:, col_ROI, 2);
    Z_aligne = ptCloud_aligne.Location(:, col_ROI, 3); 
    
    portee_min = -0.4;
    portee_max = -15.0;
    limite_gauche = -7.0;
    rho_sans_sol = rho_brut(:, col_ROI); 
    zone_lac = (X_aligne >= portee_max & X_aligne <= portee_min) & (Y_aligne >= limite_gauche);
    rho_sans_sol(~zone_lac) = 0;
    matrice_3D_ROI = cat(3, X_aligne, Y_aligne, Z_aligne);
    ptCloud_sans_sol = pointCloud(matrice_3D_ROI);
   
    % 3. SEGMENTATION DIRECTE SUR RHO BRUT
    [candidats, boites] = segmentation_rho_marin(rho_sans_sol, ptCloud_sans_sol);
    nb_objets = numel(candidats);
    
    % 4. COLORATION 
    ptCloud_visu = ptCloud_aligne; 
    [M, N, ~] = size(ptCloud_visu.Location);
    distance_max_couleur = 40.0; 
    rho_norm = rho_brut / distance_max_couleur;
    rho_norm(rho_norm > 1) = 1; 
    palette = jet(256); 
    indices_couleurs = round(rho_norm * 255) + 1;
    colors = zeros(M, N, 3, 'uint8');
    colors(:,:,1) = reshape(palette(indices_couleurs, 1) * 255, M, N);
    colors(:,:,2) = reshape(palette(indices_couleurs, 2) * 255, M, N);
    colors(:,:,3) = reshape(palette(indices_couleurs, 3) * 255, M, N);
    
    X_visu = ptCloud_visu.Location(:,:,1);
    Y_visu = ptCloud_visu.Location(:,:,2);
    Z_visu = ptCloud_visu.Location(:,:,3);
    
    if ~isempty(boites)
        for b = 1:size(boites, 1)
            box = boites(b, :);
            xMin = box(1) - box(4)/2; xMax = box(1) + box(4)/2;
            yMin = box(2) - box(5)/2; yMax = box(2) + box(5)/2;
            zMin = box(3) - box(6)/2; zMax = box(3) + box(6)/2;
            
            masque_boite = (X_visu >= xMin & X_visu <= xMax) & ...
                           (Y_visu >= yMin & Y_visu <= yMax) & ...
                           (Z_visu >= zMin & Z_visu <= zMax);
            
            canal_R = colors(:,:,1); canal_G = colors(:,:,2); canal_B = colors(:,:,3);
            canal_R(masque_boite) = 0;   canal_G(masque_boite) = 0;   canal_B(masque_boite) = 255; 
            colors(:,:,1) = canal_R; colors(:,:,2) = canal_G; colors(:,:,3) = canal_B;
        end
    end
    ptCloud_visu.Color = colors;
    
    % =====================================================================
    % 5. APPEL DU FILTRE SÉLECTIONNÉ
    % =====================================================================
    matrice_centres = [];
    if nb_objets > 0
        matrice_centres = reshape([candidats.centre], 3, [])';
    end
    
    if strcmp(CHOIX_FILTRE, 'PARTICULAIRE')
        [tracks, detections_associees, next_display_id, next_id] = filtre_particule(tracks, matrice_centres, dt, N_particules, Q_cov, R_cov, next_display_id, next_id);
  elseif strcmp(CHOIX_FILTRE, 'BOX_PARTICULAIRE')
        [tracks, detections_associees, next_display_id, next_id] = filtre_particule_boite(tracks, matrice_centres, dt, N_boites, V_modele, erreur_lidar_max, R_cov, next_display_id, next_id);
    end
        % MÉTRIQUE DE CARDINALITÉ : nombre de pistes confirmées à cette trame
    if isempty(tracks)
        m_estime = 0;
    else
        m_estime = sum([tracks.confirmed]);
    end
    historique_cardinalite = [historique_cardinalite; m_estime];
   
    % 6. AFFICHAGE 3D ET EXTRACTION
    X_final = ptCloud_visu.Location(:,:,1); Y_final = ptCloud_visu.Location(:,:,2); Z_final = ptCloud_visu.Location(:,:,3);
    idx_valides = find(X_final >= -50 & X_final <= 50 & Y_final >= -50 & Y_final <= 50 & Z_final >= -5 & Z_final <= 5);
    ptCloud_visu = select(ptCloud_visu, idx_valides); 
    
    view(viewer, ptCloud_visu);
    
    delete(findobj(viewer.Axes, 'Type', 'Patch'));
    delete(findobj(viewer.Axes, 'Tag', 'ParticulesMOT'));
    delete(findobj(viewer.Axes, 'Tag', 'LigneTrajectoireMOT'));
    delete(findobj(viewer.Axes, 'Tag', 'PointActuelMOT'));
    delete(findobj(viewer.Axes, 'Tag', 'TexteIDMOT'));
    
    if ~isempty(boites)
        showShape('cuboid', boites, 'Parent', viewer.Axes, 'Color', repmat([0 1 0], size(boites, 1), 1), 'Opacity', 0.1, 'LineWidth', 0.3);
    end
    
    if ~isempty(tracks)
        hold(viewer.Axes, 'on');
        for t = 1:numel(tracks)
            if tracks(t).confirmed
                
                % EXTRACTION DES DONNÉES SELON LE FILTRE
               % EXTRACTION DES DONNÉES SELON LE FILTRE
                if strcmp(CHOIX_FILTRE, 'PARTICULAIRE')
                    pos_actuelle = tracks(t).X_k(1:3);
                    Vx_moy = tracks(t).X_k(4); Vy_moy = tracks(t).X_k(5); Vz_moy = tracks(t).X_k(6);
                    ax_moy = tracks(t).X_k(7); ay_moy = tracks(t).X_k(8); az_moy = tracks(t).X_k(9);
                    vitesse_moy = norm([Vx_moy, Vy_moy, Vz_moy]);
                    
                    std_vx = std(tracks(t).particules(4,:)); std_vy = std(tracks(t).particules(5,:)); std_vz = std(tracks(t).particules(6,:));
                    std_ax = std(tracks(t).particules(7,:)); std_ay = std(tracks(t).particules(8,:)); std_az = std(tracks(t).particules(9,:));
                    
                    Vx_min = Vx_moy - 2*std_vx; Vx_max = Vx_moy + 2*std_vx;
                    Vy_min = Vy_moy - 2*std_vy; Vy_max = Vy_moy + 2*std_vy;
                    Vz_min = Vz_moy - 2*std_vz; Vz_max = Vz_moy + 2*std_vz;
                    
                    ax_min = ax_moy - 2*std_ax; ax_max = ax_moy + 2*std_ax;
                    ay_min = ay_moy - 2*std_ay; ay_max = ay_moy + 2*std_ay;
                    az_min = az_moy - 2*std_az; az_max = az_moy + 2*std_az;
                    
                    v_min = max(0, vitesse_moy - 2*norm([std_vx, std_vy, std_vz]));
                    v_max = vitesse_moy + 2*norm([std_vx, std_vy, std_vz]);
                    
                    texte_affichage = sprintf('ID %d | %.2f m/s', tracks(t).display_id, vitesse_moy);
                    
                    % Affichage nuage de points
                    scatter3(viewer.Axes, tracks(t).particules(1,:), tracks(t).particules(2,:), tracks(t).particules(3,:), ...
                             2, 'w', 'filled', 'MarkerFaceAlpha', 0.5, 'Tag', 'ParticulesMOT');
                             
              elseif strcmp(CHOIX_FILTRE, 'BOX_PARTICULAIRE')
                    X_hull_inf = min(inf(tracks(t).X_boxes), [], 2);
                    X_hull_sup = max(sup(tracks(t).X_boxes), [], 2);
                    X_hull = infsup(X_hull_inf, X_hull_sup);
                    
                    pos_actuelle = mid(X_hull(1:3));
                    
                    Vx_int = X_hull(4); Vy_int = X_hull(5); Vz_int = X_hull(6);
                    Ax_int = X_hull(7); Ay_int = X_hull(8); Az_int = X_hull(9);
                    V_norm_int = sqrt(Vx_int^2 + Vy_int^2 + Vz_int^2); 
                    
                    v_min = inf(V_norm_int); v_max = sup(V_norm_int); 
                    vitesse_moy = sqrt(mid(Vx_int)^2 + mid(Vy_int)^2 + mid(Vz_int)^2);
                    
                    Vx_min = inf(Vx_int); Vx_max = sup(Vx_int); Vx_moy = mid(Vx_int);
                    Vy_min = inf(Vy_int); Vy_max = sup(Vy_int); Vy_moy = mid(Vy_int);
                    Vz_min = inf(Vz_int); Vz_max = sup(Vz_int); Vz_moy = mid(Vz_int);
                    
                    ax_min = inf(Ax_int); ax_max = sup(Ax_int); ax_moy = mid(Ax_int);
                    ay_min = inf(Ay_int); ay_max = sup(Ay_int); ay_moy = mid(Ay_int);
                    az_min = inf(Az_int); az_max = sup(Az_int); az_moy = mid(Az_int);
                    
                    texte_affichage = sprintf('ID %d | V:[%.2f, %.2f]', tracks(t).display_id, v_min, v_max);
                    
                    % Affichage des Boîtes
                    for b = 1:N_boites
                        x_min = inf(tracks(t).X_boxes(1,b)); x_max = sup(tracks(t).X_boxes(1,b));
                        y_min = inf(tracks(t).X_boxes(2,b)); y_max = sup(tracks(t).X_boxes(2,b));
                        z_min = inf(tracks(t).X_boxes(3,b)); z_max = sup(tracks(t).X_boxes(3,b));
                        w = max(0.01, x_max - x_min); h = max(0.01, y_max - y_min); d = max(0.01, z_max - z_min);
                        cx = x_min + w/2; cy = y_min + h/2; cz = z_min + d/2;
                        opacite = min(1, max(0.1, tracks(t).poids(b) * 2));
                        if all(isfinite([cx, cy, cz, w, h, d]))
                            showShape('cuboid', [cx, cy, cz, w, h, d, 0, 0, 0], 'Parent', viewer.Axes, 'Color', [1 0.5 0], 'Opacity', opacite, 'LineWidth', 0.5);
                        end
                    end
                end
                
                % SAUVEGARDE ARCHIVES EN MÉMOIRE
                idx_archive = find([archive_vitesses.display_id] == tracks(t).display_id);
                bool_observe = (tracks(t).lost_frames == 0); 
                
                if isempty(idx_archive)
                    archive_vitesses = [archive_vitesses, struct('display_id', tracks(t).display_id, ...
                        'vitesses_moy', vitesse_moy, 'vitesses_min', v_min, 'vitesses_max', v_max, ...
                        'vx_moy', Vx_moy, 'vx_min', Vx_min, 'vx_max', Vx_max, ...
                        'vy_moy', Vy_moy, 'vy_min', Vy_min, 'vy_max', Vy_max, ...
                        'vz_moy', Vz_moy, 'vz_min', Vz_min, 'vz_max', Vz_max, ...
                        'ax_moy', ax_moy, 'ax_min', ax_min, 'ax_max', ax_max, ...
                        'ay_moy', ay_moy, 'ay_min', ay_min, 'ay_max', ay_max, ...
                        'az_moy', az_moy, 'az_min', az_min, 'az_max', az_max, ...
                        'frames', frame_idx, 'est_observe', bool_observe)];
                else
                    archive_vitesses(idx_archive).vitesses_moy = [archive_vitesses(idx_archive).vitesses_moy; vitesse_moy];
                    archive_vitesses(idx_archive).vitesses_min = [archive_vitesses(idx_archive).vitesses_min; v_min];
                    archive_vitesses(idx_archive).vitesses_max = [archive_vitesses(idx_archive).vitesses_max; v_max];
                    
                    archive_vitesses(idx_archive).vx_moy = [archive_vitesses(idx_archive).vx_moy; Vx_moy];
                    archive_vitesses(idx_archive).vx_min = [archive_vitesses(idx_archive).vx_min; Vx_min];
                    archive_vitesses(idx_archive).vx_max = [archive_vitesses(idx_archive).vx_max; Vx_max];
                    
                    archive_vitesses(idx_archive).vy_moy = [archive_vitesses(idx_archive).vy_moy; Vy_moy];
                    archive_vitesses(idx_archive).vy_min = [archive_vitesses(idx_archive).vy_min; Vy_min];
                    archive_vitesses(idx_archive).vy_max = [archive_vitesses(idx_archive).vy_max; Vy_max];
                    
                    archive_vitesses(idx_archive).vz_moy = [archive_vitesses(idx_archive).vz_moy; Vz_moy];
                    archive_vitesses(idx_archive).vz_min = [archive_vitesses(idx_archive).vz_min; Vz_min];
                    archive_vitesses(idx_archive).vz_max = [archive_vitesses(idx_archive).vz_max; Vz_max];
                    
                    archive_vitesses(idx_archive).ax_moy = [archive_vitesses(idx_archive).ax_moy; ax_moy];
                    archive_vitesses(idx_archive).ax_min = [archive_vitesses(idx_archive).ax_min; ax_min];
                    archive_vitesses(idx_archive).ax_max = [archive_vitesses(idx_archive).ax_max; ax_max];
                    
                    archive_vitesses(idx_archive).ay_moy = [archive_vitesses(idx_archive).ay_moy; ay_moy];
                    archive_vitesses(idx_archive).ay_min = [archive_vitesses(idx_archive).ay_min; ay_min];
                    archive_vitesses(idx_archive).ay_max = [archive_vitesses(idx_archive).ay_max; ay_max];
                    
                    archive_vitesses(idx_archive).az_moy = [archive_vitesses(idx_archive).az_moy; az_moy];
                    archive_vitesses(idx_archive).az_min = [archive_vitesses(idx_archive).az_min; az_min];
                    archive_vitesses(idx_archive).az_max = [archive_vitesses(idx_archive).az_max; az_max];

                    archive_vitesses(idx_archive).frames = [archive_vitesses(idx_archive).frames; frame_idx];
                    archive_vitesses(idx_archive).est_observe = [archive_vitesses(idx_archive).est_observe; bool_observe];
                end
                
                % DESSIN TRAJECTOIRES
                if size(tracks(t).historique, 1) > 1
                    if strcmp(CHOIX_FILTRE, 'PARTICULAIRE')
                        plot3(viewer.Axes, tracks(t).historique(:,1), tracks(t).historique(:,2), tracks(t).historique(:,3), 'c-', 'LineWidth', 2.0, 'Tag', 'LigneTrajectoireMOT');
                    else
                        plot3(viewer.Axes, tracks(t).historique(:,1), tracks(t).historique(:,2), tracks(t).historique(:,3), '-', 'Color', [1 0.5 0], 'LineWidth', 2.0, 'Tag', 'LigneTrajectoireMOT');
                    end
                end
                
                % Le point central 
                if strcmp(CHOIX_FILTRE, 'PARTICULAIRE')
                    plot3(viewer.Axes, pos_actuelle(1), pos_actuelle(2), pos_actuelle(3), 'co', 'MarkerSize', 5.0, 'MarkerFaceColor', 'c', 'Tag', 'PointActuelMOT');
                else
                    plot3(viewer.Axes, pos_actuelle(1), pos_actuelle(2), pos_actuelle(3), 'o', 'Color', [1 0.5 0], 'MarkerSize', 4.0, 'MarkerFaceColor', [1 0.5 0], 'Tag', 'PointActuelMOT');
                end
                text(viewer.Axes, pos_actuelle(1), pos_actuelle(2), pos_actuelle(3) + 0.8, texte_affichage, 'Color', 'white', 'FontSize', 8, 'FontWeight', 'bold', 'Tag', 'TexteIDMOT');    
            end 
       end
        hold(viewer.Axes, 'off');   
    end
    drawnow limitrate;
end
disp('Fin de l''acquisition marine.');
% =========================================================================
% MÉTRIQUE D'ÉVALUATION MOT : ERREUR DE CARDINALITÉ
% =========================================================================
if ~isempty(historique_cardinalite)
    frames_card = (1:numel(historique_cardinalite))';
    erreur_cardinalite = abs(n_reel - historique_cardinalite);
    
    fprintf('\n--- ERREUR DE CARDINALITE (%s) ---\n', CHOIX_FILTRE);
    fprintf('Nombre reel de cibles          : %d\n', n_reel);
    fprintf('Nombre estime moyen            : %.2f\n', mean(historique_cardinalite));
    fprintf('Erreur de cardinalite moyenne  : %.2f\n', mean(erreur_cardinalite));
    fprintf('Trames a cardinalite exacte    : %.1f %%\n', 100 * sum(erreur_cardinalite == 0) / numel(erreur_cardinalite));
    fprintf('Erreur de cardinalite maximale : %d\n', max(erreur_cardinalite));
    disp('--------------------------------------------------');
    
    figure('Name', ['Erreur de Cardinalite - ' CHOIX_FILTRE], 'NumberTitle', 'off', 'Color', 'w');
    
    subplot(2,1,1); hold on;
    plot(frames_card, n_reel*ones(size(frames_card)), '--k', 'LineWidth', 1.5, 'DisplayName', 'Nombre reel de cibles');
    plot(frames_card, historique_cardinalite, '-b', 'LineWidth', 1.5, 'DisplayName', 'Nombre estime de cibles');
    hold off; grid on;
    ylabel('Nombre de cibles', 'FontWeight', 'bold');
    title(['Estimation de la cardinalite - ' CHOIX_FILTRE], 'FontSize', 12);
    legend('Location', 'best');
    ylim([0, max(n_reel, max(historique_cardinalite)) + 1]);
    
    subplot(2,1,2);
    plot(frames_card, erreur_cardinalite, '-r', 'LineWidth', 1.5);
    grid on;
    xlabel('Temps (Frames)', 'FontWeight', 'bold');
    ylabel('Erreur |n - m|', 'FontWeight', 'bold');
    title('Erreur de cardinalite au cours du temps', 'FontSize', 12);
    ylim([0, max(erreur_cardinalite) + 1]);
end
% =====================================================
% 7. GÉNÉRATION DES GRAPHIQUES DE VITESSE
% =====================================================
if ~isempty(archive_vitesses)
    
    archive_vitesses_filtree = [];
    for k = 1:numel(archive_vitesses)
        duree_vie = length(archive_vitesses(k).frames);
        nb_observations_reelles = sum(archive_vitesses(k).est_observe);
        ratio_observation = nb_observations_reelles / duree_vie;
        
        if duree_vie > 10 && ratio_observation > 0.9
            fprintf('✅ ID %d validé : Durée = %d frames | Ratio = %.2f\n', archive_vitesses(k).display_id, duree_vie, ratio_observation);
            archive_vitesses_filtree = [archive_vitesses_filtree, archive_vitesses(k)];
        else
            fprintf('❌ ID %d rejeté (Fantôme) : Durée = %d frames | Ratio = %.2f\n', archive_vitesses(k).display_id, duree_vie, ratio_observation);
        end
    end
   
    % GRAPHIQUE 1 : LA VITESSE GLOBALE (AVANT supression des anomalie)
    figure('Name', sprintf('Avant Supression des anomalies - %s', CHOIX_FILTRE), 'NumberTitle', 'off', 'Color', 'w');
    hold on;
    couleurs_plot_brut = hsv(numel(archive_vitesses)); 
    
    for k = 1:numel(archive_vitesses)
        v = archive_vitesses(k).vitesses_moy; 
        f = archive_vitesses(k).frames;
        plot(f, movmean(v, fenetre_lissage), '-o', 'MarkerSize', 6, 'Color', couleurs_plot_brut(k,:), 'DisplayName', sprintf('ID %d', archive_vitesses(k).display_id));
    end
    hold off; grid on;
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('Vitesse Absolue (m/s)', 'FontWeight', 'bold');
    title(sprintf('Évolution de la vitesse globale (Avec Fantômes) - %s', CHOIX_FILTRE), 'FontSize', 12);
    legend('Location', 'eastoutside', 'FontSize', 8, 'NumColumns', 2);
    
    % GRAPHIQUE 2 : LA VITESSE GLOBALE (APRÈS Supression des anomalies )
    if ~isempty(archive_vitesses_filtree)
        figure('Name', sprintf('Après Supression des anomalies - %s', CHOIX_FILTRE), 'NumberTitle', 'off', 'Color', 'w');
        hold on;
        couleurs_plot_net = lines(numel(archive_vitesses_filtree)); 
        
        for k = 1:numel(archive_vitesses_filtree)
            v_moy = archive_vitesses_filtree(k).vitesses_moy;
            v_min = max(0, archive_vitesses_filtree(k).vitesses_min); 
            v_max = archive_vitesses_filtree(k).vitesses_max;
            f = archive_vitesses_filtree(k).frames;
            
            moyenne_globale = mean(v_moy); 
            
            X_fill = [f', fliplr(f')];
            Y_fill = [v_min', fliplr(v_max')];
            % Affichage du couloir d'incertitude
            fill(X_fill, Y_fill, couleurs_plot_net(k,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(f, movmean(v_moy, fenetre_lissage), '-', 'LineWidth', 1.5, 'Color', couleurs_plot_net(k,:), 'DisplayName', sprintf('ID %d (Moy: %.2f m/s)', archive_vitesses_filtree(k).display_id, moyenne_globale));
            plot([min(f), max(f)], [moyenne_globale, moyenne_globale], '--', 'LineWidth', 2, 'Color', couleurs_plot_net(k,:), 'HandleVisibility', 'off');
        end
        hold off; grid on;
        xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('Vitesse Absolue (m/s)', 'FontWeight', 'bold');
        title(sprintf('Évolution de la vitesse globale (%s)', CHOIX_FILTRE), 'FontSize', 12);
        legend('Location', 'best', 'FontSize', 10);
        
        % GRAPHIQUE 3 : DÉCOMPOSITION (Vx, Vy, Vz)
        for k = 1:numel(archive_vitesses_filtree)
            id_obj = archive_vitesses_filtree(k).display_id;
            f  = archive_vitesses_filtree(k).frames;
            
            vx_moy = archive_vitesses_filtree(k).vx_moy;
            vy_moy = archive_vitesses_filtree(k).vy_moy;
            vz_moy = archive_vitesses_filtree(k).vz_moy;
            
            figure('Name', sprintf('Décomposition Vitesse - ID %d', id_obj), 'NumberTitle', 'off', 'Color', 'w');
            hold on;
            
            % DESSIN DES COULOIRS D'INCERTITUDE
            X_fill = [f', fliplr(f')];
            
            Y_fill_x = [archive_vitesses_filtree(k).vx_min', fliplr(archive_vitesses_filtree(k).vx_max')];
            fill(X_fill, Y_fill_x, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            
            Y_fill_y = [archive_vitesses_filtree(k).vy_min', fliplr(archive_vitesses_filtree(k).vy_max')];
            fill(X_fill, Y_fill_y, 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            
            Y_fill_z = [archive_vitesses_filtree(k).vz_min', fliplr(archive_vitesses_filtree(k).vz_max')];
            fill(X_fill, Y_fill_z, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            
            plot(f, movmean(vx_moy, fenetre_lissage), '-r', 'LineWidth', 1.5, 'DisplayName', 'Vx (Axe d''approche X)');
            plot(f, movmean(vy_moy, fenetre_lissage), '-g', 'LineWidth', 1.5, 'DisplayName', 'Vy (Axe latéral Y)');
            plot(f, movmean(vz_moy, fenetre_lissage), '-b', 'LineWidth', 1.5, 'DisplayName', 'Vz (Axe vertical Z)');
            plot([min(f), max(f)], [0, 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            hold off; grid on;
            xlabel('Temps (Frames)', 'FontWeight', 'bold'); 
            ylabel('Vitesse (m/s)', 'FontWeight', 'bold');
            title(sprintf('Composantes 3D de la vitesse (Vx, Vy, Vz) - ID %d', id_obj), 'FontSize', 12);
            legend('Location', 'best', 'FontSize', 10);
        end
      % GRAPHIQUE 4 : DÉCOMPOSITION ACCÉLÉRATION (ax, ay, az)
        for k = 1:numel(archive_vitesses_filtree)
            id_obj = archive_vitesses_filtree(k).display_id;
            f  = archive_vitesses_filtree(k).frames;
            ax_moy = archive_vitesses_filtree(k).ax_moy; ay_moy = archive_vitesses_filtree(k).ay_moy; az_moy = archive_vitesses_filtree(k).az_moy;
            
            figure('Name', sprintf('Décomposition Accélération - ID %d', id_obj), 'NumberTitle', 'off', 'Color', 'w');
            hold on;
            
            X_fill = [f', fliplr(f')];
            Y_fill_x = [archive_vitesses_filtree(k).ax_min', fliplr(archive_vitesses_filtree(k).ax_max')];
            fill(X_fill, Y_fill_x, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            Y_fill_y = [archive_vitesses_filtree(k).ay_min', fliplr(archive_vitesses_filtree(k).ay_max')];
            fill(X_fill, Y_fill_y, 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            Y_fill_z = [archive_vitesses_filtree(k).az_min', fliplr(archive_vitesses_filtree(k).az_max')];
            fill(X_fill, Y_fill_z, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            
            plot(f, movmean(ax_moy, fenetre_lissage), '-r', 'LineWidth', 1.5, 'DisplayName', 'ax (Axe X)');
            plot(f, movmean(ay_moy, fenetre_lissage), '-g', 'LineWidth', 1.5, 'DisplayName', 'ay (Axe Y)');
            plot(f, movmean(az_moy, fenetre_lissage), '-b', 'LineWidth', 1.5, 'DisplayName', 'az (Axe Z)');
            plot([min(f), max(f)], [0, 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            hold off; grid on;
            xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('Accélération (m/s²)', 'FontWeight', 'bold');
            title(sprintf('Accélérations Estimées (ax, ay, az) - OBJET ID %d', id_obj), 'FontSize', 12);
            legend('Location', 'best', 'FontSize', 10);
        end  
    else
        disp('Aucun objet n''a passé le filtre de ratio > 0.5. (Graphiques propres non générés)');
    end
else
    disp('Aucun objet n''a été détecté durant la session.');
end