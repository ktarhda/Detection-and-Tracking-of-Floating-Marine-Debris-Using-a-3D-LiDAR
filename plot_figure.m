clear; clc; close all;
disp('  ÉTUDE COMPARATIVE : CLASSIQUE vs ENSEMBLISTE vs PARTICULAIRE vs BOX PARTICULAIRE  ');
fenetre_lissage = 9;   % Taille de la fenêtre du filtre moyenne glissante (en frames)
% =========================================================================
% 1. CHARGEMENT DES 4 FICHIERS .MAT
% =========================================================================
try
    data_C = load('donnees_CLASSIQUE.mat');
    archive_C = data_C.archive_vitesses_filtree;
    disp('✅ Fichier Classique chargé.');
catch
    error('Fichier donnees_CLASSIQUE.mat introuvable. Lancez le sauvegarde.m avec le filtre CLASSIQUE d''abord.');
end

try
    data_E = load('donnees_ENSEMBLISTE.mat');
    archive_E = data_E.archive_vitesses_filtree;
    disp('✅ Fichier Ensembliste chargé.');
catch
    error('Fichier donnees_ENSEMBLISTE.mat introuvable. Lancez le sauvegarde.m avec le filtre ENSEMBLISTE d''abord.');
end

try
    data_P = load('donnees_PARTICULAIRE.mat');
    archive_P = data_P.archive_vitesses_filtree;
    disp('✅ Fichier Particulaire chargé.');
catch
    error('Fichier donnees_PARTICULAIRE.mat introuvable. Lancez le sauvegarde.m avec le filtre PARTICULAIRE d''abord.');
end

try
    data_BP = load('donnees_BOX_PARTICULAIRE.mat');
    archive_BP = data_BP.archive_vitesses_filtree;
    disp('✅ Fichier Box Particulaire chargé.');
catch
    error('Fichier donnees_BOX_PARTICULAIRE.mat introuvable. Lancez le sauvegarde.m avec le filtre BOX_PARTICULAIRE d''abord.');
end

% =========================================================================
% 2. RECHERCHE DES OBJETS COMMUNS (Aux 4 méthodes)
% =========================================================================
IDs_C = [archive_C.display_id];
IDs_E = [archive_E.display_id];
IDs_P = [archive_P.display_id];
IDs_BP = [archive_BP.display_id];

% On trouve les ID qui existent dans les 4 méthodes
IDs_communs = intersect(intersect(intersect(IDs_C, IDs_E), IDs_P), IDs_BP);

if isempty(IDs_communs)
    disp('❌ Aucun ID d''objet commun trouvé entre les quatre exécutions.');
    return;
else
    fprintf('🎯 %d objet(s) commun(s) trouvé(s) pour la comparaison !\n', length(IDs_communs));
end

% Couleurs pour la comparaison
c_classique = 'r';            % Rouge pour Classique
c_ensembliste = 'g';          % Vert pour Ensembliste
c_particulaire = 'b';         % Bleu pour Particulaire
c_box_particulaire = [1 0.5 0]; % Orange pour Box Particulaire

% =========================================================================
% 3. GÉNÉRATION DES GRAPHIQUES SUPERPOSÉS
% =========================================================================
for i = 1:length(IDs_communs)
    id_obj = IDs_communs(i);
    
    % Extraction des données de l'objet pour les 4 méthodes
    obj_C = archive_C(IDs_C == id_obj);
    obj_E = archive_E(IDs_E == id_obj);
    obj_P = archive_P(IDs_P == id_obj);
    obj_BP = archive_BP(IDs_BP == id_obj);
    
    % --- SYNCHRONISATION TEMPORELLE STRICTE (4 Filtres) ---
    [frames_temp1, iC_t1, iE_t1] = intersect(obj_C.frames, obj_E.frames);
    [frames_temp2, iCE_t2, iP_t2] = intersect(frames_temp1, obj_P.frames);
    [frames_communes, iCEP_t3, idx_BP] = intersect(frames_temp2, obj_BP.frames);
    
    idx_C = iC_t1(iCE_t2(iCEP_t3));
    idx_E = iE_t1(iCE_t2(iCEP_t3));
    idx_P = iP_t2(iCEP_t3);
    
    if isempty(frames_communes)
        continue; 
    end

    %% GRAPHIQUE 1 : TRAJECTOIRE 3D SUPERPOSÉE
    figure('Name', sprintf('Comparaison Trajectoire 3D - ID %d', id_obj), 'Color', 'w');
    hold on;
    plot3(movmean(obj_C.x(idx_C), fenetre_lissage), movmean(obj_C.y(idx_C), fenetre_lissage), movmean(obj_C.z(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Filtre Classique');
    plot3(movmean(obj_E.x(idx_E), fenetre_lissage), movmean(obj_E.y(idx_E), fenetre_lissage), movmean(obj_E.z(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0, 'DisplayName', 'Filtre Ensembliste');
    plot3(movmean(obj_P.x(idx_P), fenetre_lissage), movmean(obj_P.y(idx_P), fenetre_lissage), movmean(obj_P.z(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Filtre Particulaire');
    plot3(movmean(obj_BP.x(idx_BP), fenetre_lissage), movmean(obj_BP.y(idx_BP), fenetre_lissage), movmean(obj_BP.z(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0, 'DisplayName', 'Box Particulaire');
    
    % Point de départ
    plot3(obj_C.x(idx_C(1)), obj_C.y(idx_C(1)), obj_C.z(idx_C(1)), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Départ');
    
    grid on; view(3);
    xlabel('Position X (m)', 'FontWeight', 'bold'); ylabel('Position Y (m)', 'FontWeight', 'bold'); zlabel('Position Z (m)', 'FontWeight', 'bold');
    title(sprintf('Comparaison des Trajectoires 3D - Objet ID %d', id_obj), 'FontSize', 12);
    legend('Location', 'best');
    
        %% GRAPHIQUE 2 : VITESSE GLOBALE
    figure('Name', sprintf('Comparaison Vitesse Globale - ID %d', id_obj), 'Color', 'w');
    hold on;
    plot(frames_communes, movmean(obj_C.vitesses_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Filtre Classique');
    plot(frames_communes, movmean(obj_E.vitesses_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0, 'DisplayName', 'Filtre Ensembliste');
    plot(frames_communes, movmean(obj_P.vitesses_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Filtre Particulaire');
    plot(frames_communes, movmean(obj_BP.vitesses_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0, 'DisplayName', 'Box Particulaire');
    
    grid on;
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('Vitesse Absolue (m/s)', 'FontWeight', 'bold');
    title(sprintf('Comparaison de la Vitesse Globale Estimée - Objet ID %d', id_obj), 'FontSize', 12);
    legend('Location', 'best');
    
      %% GRAPHIQUE 3 : COMPOSANTES DE LA VITESSE (Vx, Vy, Vz)
    figure('Name', sprintf('Vitesses Vx, Vy, Vz - ID %d', id_obj), 'Color', 'w');
    
    subplot(3,1,1); hold on;
    plot(frames_communes, movmean(obj_C.vx_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Classique');
    plot(frames_communes, movmean(obj_E.vx_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 1.5, 'DisplayName', 'Ensembliste');
    plot(frames_communes, movmean(obj_P.vx_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Particulaire');
    plot(frames_communes, movmean(obj_BP.vx_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Box Particulaire');
    ylabel('Vx (m/s)', 'FontWeight', 'bold'); title(sprintf('Décomposition Vitesse - ID %d', id_obj)); grid on; legend('Location', 'best');
    
    subplot(3,1,2); hold on;
    plot(frames_communes, movmean(obj_C.vy_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.vy_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_P.vy_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.vy_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 1.5);
    ylabel('Vy (m/s)', 'FontWeight', 'bold'); grid on;
    
    subplot(3,1,3); hold on;
    plot(frames_communes, movmean(obj_C.vz_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.vz_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_P.vz_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.vz_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 1.5);
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('Vz (m/s)', 'FontWeight', 'bold'); grid on;
    
     %% GRAPHIQUE 4 : ÉVOLUTION DE L'INCERTITUDE DE POSITION (P11, P22, P33)
    figure('Name', sprintf('Incertitude Position - ID %d', id_obj), 'Color', 'w');
    
    subplot(3,1,1); hold on;
    plot(frames_communes, movmean(obj_C.p11(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Classique P_{11}');
    plot(frames_communes, movmean(obj_E.p11(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0, 'DisplayName', 'Ensembliste P^+_{11}');
    plot(frames_communes, movmean(obj_P.p11(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Particulaire Var(X)');
    plot(frames_communes, movmean(obj_BP.p11(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0, 'DisplayName', 'Box Particulaire (Largeur)');
    ylabel('P_{11} (X)', 'FontWeight', 'bold'); 
    title(sprintf('Comparaison des Covariances de Position - ID %d', id_obj)); grid on; legend('Location', 'best');
    
    subplot(3,1,2); hold on;
    plot(frames_communes, movmean(obj_C.p22(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.p22(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0);
    plot(frames_communes, movmean(obj_P.p22(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.p22(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0);
    ylabel('P_{22} (Y)', 'FontWeight', 'bold'); grid on;
    
    subplot(3,1,3); hold on;
    plot(frames_communes, movmean(obj_C.p33(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.p33(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0);
    plot(frames_communes, movmean(obj_P.p33(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.p33(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0);
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('P_{33} (Z)', 'FontWeight', 'bold'); grid on;
      %% GRAPHIQUE 5 : ÉVOLUTION DE L'INCERTITUDE DE VITESSE (P44, P55, P66)
    figure('Name', sprintf('Incertitude Vitesse - ID %d', id_obj), 'Color', 'w');
    
    subplot(3,1,1); hold on;
    plot(frames_communes, movmean(obj_C.p44(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Classique P_{44}');
    plot(frames_communes, movmean(obj_E.p44(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0, 'DisplayName', 'Ensembliste P^+_{44}');
    plot(frames_communes, movmean(obj_P.p44(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Particulaire Var(Vx)');
    plot(frames_communes, movmean(obj_BP.p44(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0, 'DisplayName', 'Box Particulaire (Largeur)');
    ylabel('P_{44} (Vx)', 'FontWeight', 'bold'); 
    title(sprintf('Comparaison des Covariances de Vitesse - ID %d', id_obj)); grid on; legend('Location', 'best');
    
    subplot(3,1,2); hold on;
    plot(frames_communes, movmean(obj_C.p55(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.p55(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0);
    plot(frames_communes, movmean(obj_P.p55(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.p55(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0);
    ylabel('P_{55} (Vy)', 'FontWeight', 'bold'); grid on;
    
    subplot(3,1,3); hold on;
    plot(frames_communes, movmean(obj_C.p66(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.p66(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0);
    plot(frames_communes, movmean(obj_P.p66(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.p66(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0);
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('P_{66} (Vz)', 'FontWeight', 'bold'); grid on;
    
     %% GRAPHIQUE 6 : COMPOSANTES DE L'ACCÉLÉRATION (ax, ay, az)
    figure('Name', sprintf('Accélération ax, ay, az - ID %d', id_obj), 'Color', 'w');
    
    subplot(3,1,1); hold on;
    plot(frames_communes, movmean(obj_C.ax_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Classique');
    plot(frames_communes, movmean(obj_E.ax_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 1.5, 'DisplayName', 'Ensembliste');
    plot(frames_communes, movmean(obj_P.ax_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Particulaire');
    plot(frames_communes, movmean(obj_BP.ax_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Box Particulaire');
    ylabel('ax (m/s^2)', 'FontWeight', 'bold'); title(sprintf('Décomposition Accélération - ID %d', id_obj)); grid on; legend('Location', 'best');
    
    subplot(3,1,2); hold on;
    plot(frames_communes, movmean(obj_C.ay_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.ay_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_P.ay_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.ay_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 1.5);
    ylabel('ay (m/s^2)', 'FontWeight', 'bold'); grid on;
    
    subplot(3,1,3); hold on;
    plot(frames_communes, movmean(obj_C.az_moy(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.az_moy(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_P.az_moy(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.az_moy(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 1.5);
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('az (m/s^2)', 'FontWeight', 'bold'); grid on;

    %% GRAPHIQUE 7 : ÉVOLUTION DE L'INCERTITUDE D'ACCÉLÉRATION (P77, P88, P99)
    figure('Name', sprintf('Incertitude Accélération - ID %d', id_obj), 'Color', 'w');
    
    subplot(3,1,1); hold on;
    plot(frames_communes, movmean(obj_C.p77(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5, 'DisplayName', 'Classique P_{77}');
    plot(frames_communes, movmean(obj_E.p77(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0, 'DisplayName', 'Ensembliste P^+_{77}');
    plot(frames_communes, movmean(obj_P.p77(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5, 'DisplayName', 'Particulaire Var(ax)');
    plot(frames_communes, movmean(obj_BP.p77(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0, 'DisplayName', 'Box Particulaire (Largeur)');
    ylabel('P_{77} (ax)', 'FontWeight', 'bold'); 
    title(sprintf('Comparaison des Covariances d''Accélération - ID %d', id_obj)); grid on; legend('Location', 'best');
    
    subplot(3,1,2); hold on;
    plot(frames_communes, movmean(obj_C.p88(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.p88(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0);
    plot(frames_communes, movmean(obj_P.p88(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.p88(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0);
    ylabel('P_{88} (ay)', 'FontWeight', 'bold'); grid on;
    
    subplot(3,1,3); hold on;
    plot(frames_communes, movmean(obj_C.p99(idx_C), fenetre_lissage), '-', 'Color', c_classique, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_E.p99(idx_E), fenetre_lissage), '-', 'Color', c_ensembliste, 'LineWidth', 2.0);
    plot(frames_communes, movmean(obj_P.p99(idx_P), fenetre_lissage), '-', 'Color', c_particulaire, 'LineWidth', 1.5);
    plot(frames_communes, movmean(obj_BP.p99(idx_BP), fenetre_lissage), '-', 'Color', c_box_particulaire, 'LineWidth', 2.0);
    xlabel('Temps (Frames)', 'FontWeight', 'bold'); ylabel('P_{99} (az)', 'FontWeight', 'bold'); grid on;
end

disp('✅ Graphiques comparatifs des 4 filtres générés avec succès !');

% =========================================================================
% 4. MÉTRIQUE D'ÉVALUATION MOT : ERREUR DE CARDINALITÉ
% =========================================================================
noms_filtres  = {'Classique', 'Ensembliste', 'Particulaire', 'Box Particulaire'};
data_filtres  = {data_C, data_E, data_P, data_BP};
couleurs_card = {c_classique, c_ensembliste, c_particulaire, c_box_particulaire};

figure('Name', 'Comparaison Erreur de Cardinalite', 'Color', 'w');

fprintf('\n================== ERREUR DE CARDINALITE ==================\n');
fprintf('%-18s | %8s | %8s | %9s | %5s\n', 'Filtre', 'Moy est', 'Err moy', 'Exact (%)', 'Max');
fprintf('-----------------------------------------------------------\n');

n_reel_ref = NaN;
for f = 1:4
    if ~isfield(data_filtres{f}, 'historique_cardinalite')
        fprintf('%-18s | (donnees absentes : relancer sauvegarde.m)\n', noms_filtres{f});
        continue;
    end
    
    card = data_filtres{f}.historique_cardinalite;
    n_reel_ref = data_filtres{f}.n_reel;
    err = abs(n_reel_ref - card);
    fr = (1:numel(card))';
    
    fprintf('%-18s | %8.2f | %8.2f | %9.1f | %5d\n', noms_filtres{f}, ...
            mean(card), mean(err), 100*sum(err == 0)/numel(err), max(err));
    
    subplot(2,1,1); hold on;
    plot(fr, card, '-', 'Color', couleurs_card{f}, 'LineWidth', 1.5, 'DisplayName', noms_filtres{f});
    
    subplot(2,1,2); hold on;
    plot(fr, err, '-', 'Color', couleurs_card{f}, 'LineWidth', 1.5, 'DisplayName', noms_filtres{f});
end
fprintf('===========================================================\n\n');

if ~isnan(n_reel_ref)
    subplot(2,1,1);
    plot(xlim, [n_reel_ref n_reel_ref], '--k', 'LineWidth', 1.5, 'DisplayName', 'Nombre reel de cibles');
end

subplot(2,1,1);
grid on; ylabel('Nombre de cibles', 'FontWeight', 'bold');
title('Estimation de la cardinalite - Comparaison des 4 filtres', 'FontSize', 12);
legend('Location', 'best');

subplot(2,1,2);
grid on; xlabel('Temps (Frames)', 'FontWeight', 'bold');
ylabel('Erreur |n - m|', 'FontWeight', 'bold');
title('Erreur de cardinalite au cours du temps', 'FontSize', 12);
legend('Location', 'best');