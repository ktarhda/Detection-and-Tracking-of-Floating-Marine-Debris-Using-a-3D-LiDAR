"""Filtre de Kalman ensembliste UBIKF (kalman_ensembliste.m).

L'état est une boîte (vecteur de 9 intervalles), la covariance P+ est une
borne supérieure réelle (Tran et al., 2021). La correction suit Xiong et al.
(2013) : gain de Kalman, puis propagation de contraintes (intersection avec
la boîte de mesure, contraintes backward sur la vitesse et l'accélération).

Champs d'une piste : id, X_k_intval, P_k_plus, X_prec_intval, historique,
historique_observe, lost_frames, age, confirmed, display_id.
"""

import numpy as np

from ..intervalles import Intervalle, intersect, produit
from ..parametres import matrice_Q_jerk
from .commun import ajouter_ligne, confirmer, distances_mahalanobis, nb_detections


def kalman_ensembliste(tracks, matrice_centres, dt, A, C, Q, R, V_modele, erreur_lidar_max,
                       next_display_id, next_id, seuil_mahalanobis=3.5, age_confirmation=10,
                       trames_perdues_max=10, q_bruit=0.01, facteur_sigma=1.5, rayons_initiaux_va=None):
    nb_objets = nb_detections(matrice_centres)
    detections_associees = np.zeros(nb_objets, dtype=bool)
    if rayons_initiaux_va is None:
        rayons_initiaux_va = np.array([0.5, 0.5, 0.5, 0.1, 0.1, 0.1])

    # OPTIMISATION 1 : matrice Q cinématique (bruit de jerk). Q passé en argument
    # n'est pas utilisé pour P+, exactement comme dans le code MATLAB.
    Q_opt = matrice_Q_jerk(dt, q_bruit)
    I9 = np.eye(9)

    # 1. PRÉDICTION ENSEMBLISTE
    for piste in tracks:
        piste["X_prec_intval"] = piste["X_k_intval"].copy()          # pour le forward-backward
        piste["X_k_intval"] = produit(A, piste["X_k_intval"]) + V_modele
        piste["P_k_plus"] = A @ piste["P_k_plus"] @ A.T + Q_opt

    # 2. ASSOCIATION (MAHALANOBIS) ET CORRECTION
    for piste in tracks:
        if nb_objets == 0:
            break

        X = piste["X_k_intval"]
        centre_predit = X[0:3].mid()
        S_k = C @ piste["P_k_plus"] @ C.T + R
        distances = distances_mahalanobis(matrice_centres, centre_predit, np.linalg.inv(S_k),
                                          detections_associees)
        idx_best = int(np.argmin(distances))
        min_mahal = distances[idx_best]

        if min_mahal < seuil_mahalanobis:
            mesure_Y = matrice_centres[idx_best]
            boite_lidar = Intervalle.midrad(mesure_Y, erreur_lidar_max)

            # Gain optimal à partir de la borne supérieure P+
            K = np.linalg.solve(S_k.T, (piste["P_k_plus"] @ C.T).T).T
            # Équations a1, a2, a3 (Xiong et al., 2013)
            a1 = produit(C, X)
            a2 = boite_lidar - a1
            a3 = produit(K, a2)
            X_temp = X + a3
            # Mise à jour de la borne supérieure P+ (Tran et al., 2021)
            piste["P_k_plus"] = (I9 - K @ C) @ piste["P_k_plus"]

            # OPTIMISATION 2 : contrainte spatiale (intersection avec la mesure)
            nouvelle_position = intersect(X_temp[0:3], boite_lidar)
            if nouvelle_position.contient_nan():
                nouvelle_position = boite_lidar
            X[0:3] = nouvelle_position

            # Contrainte cinématique backward sur la vitesse
            V_backward = (nouvelle_position - piste["X_prec_intval"][0:3]) / dt
            sigma_v = np.sqrt(np.diag(piste["P_k_plus"])[3:6])
            V_ubikf = Intervalle.midrad(X_temp[3:6].mid(), facteur_sigma * sigma_v)
            vitesse_contractee = intersect(V_ubikf, V_backward)
            X[3:6] = V_ubikf if vitesse_contractee.contient_nan() else vitesse_contractee

            # Contrainte backward sur l'accélération
            A_backward = (X[3:6] - piste["X_prec_intval"][3:6]) / dt
            sigma_a = np.sqrt(np.diag(piste["P_k_plus"])[6:9])
            A_ubikf = Intervalle.midrad(X_temp[6:9].mid(), facteur_sigma * sigma_a)
            accel_contractee = intersect(A_ubikf, A_backward)
            X[6:9] = A_ubikf if accel_contractee.contient_nan() else accel_contractee

            piste["historique"] = ajouter_ligne(piste["historique"], X[0:3].mid())
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], mesure_Y)
            piste["lost_frames"] = 0
            piste["age"] += 1
            next_display_id = confirmer(piste, next_display_id, age_confirmation)
            detections_associees[idx_best] = True
        else:
            # COASTING : vitesse et accélération gelées, avec une largeur de 1,5 sigma
            sigma_v = np.sqrt(np.diag(piste["P_k_plus"])[3:6])
            X[3:6] = Intervalle.midrad(X[3:6].mid(), facteur_sigma * sigma_v)
            sigma_a = np.sqrt(np.diag(piste["P_k_plus"])[6:9])
            X[6:9] = Intervalle.midrad(X[6:9].mid(), facteur_sigma * sigma_a)
            piste["lost_frames"] += 1
            piste["historique"] = ajouter_ligne(piste["historique"], X[0:3].mid())
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], [np.nan] * 3)

    # 3. NETTOYAGE
    tracks = [p for p in tracks if p["lost_frames"] <= trames_perdues_max]

    # 4. CRÉATION D'UNE NOUVELLE PISTE
    for d in range(nb_objets):
        if not detections_associees[d]:
            centre = matrice_centres[d]
            centre_initial = np.concatenate([centre, np.zeros(6)])
            incertitude = np.concatenate([np.full(3, erreur_lidar_max), rayons_initiaux_va])
            tracks.append({
                "id": next_id,
                "X_k_intval": Intervalle.midrad(centre_initial, incertitude),
                "P_k_plus": np.eye(9),
                "X_prec_intval": Intervalle.midrad(centre_initial, incertitude),
                "historique": centre.reshape(1, 3).copy(),
                "historique_observe": centre.reshape(1, 3).copy(),
                "lost_frames": 0,
                "age": 1,
                "confirmed": False,
                "display_id": 0,
            })
            next_id += 1

    return tracks, detections_associees, next_display_id, next_id
