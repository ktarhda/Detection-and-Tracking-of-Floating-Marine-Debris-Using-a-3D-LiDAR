"""Filtre de Kalman classique multi-objets (kalman_classique.m).

Une piste est un dict avec les mêmes champs que la structure MATLAB :
id, X_k (9,), P_k (9, 9), historique (n, 3), historique_observe (n, 3),
lost_frames, age, confirmed, display_id.
"""

import numpy as np

from .commun import ajouter_ligne, confirmer, nb_detections, sans_nan


def kalman_classique(tracks, matrice_centres, dt, A, C, Q, R, next_display_id, next_id,
                     seuil_mahalanobis=3.5, age_confirmation=5, trames_perdues_max=10, P_initiale=10.0):
    nb_objets = nb_detections(matrice_centres)
    detections_associees = np.zeros(nb_objets, dtype=bool)

    # 1. PRÉDICTION
    for piste in tracks:
        piste["X_k"] = A @ piste["X_k"]
        piste["P_k"] = A @ piste["P_k"] @ A.T + Q

    # 2. ASSOCIATION ET CORRECTION (Mahalanobis)
    for piste in tracks:
        if nb_objets == 0:
            break  # comme MATLAB : sans détection, les pistes ne sont pas comptées comme perdues

        S = C @ piste["P_k"] @ C.T + R
        invS = np.linalg.inv(S)
        diff = matrice_centres - piste["X_k"][:3]
        distances = np.sqrt(np.sum((diff @ invS) * diff, axis=1))
        distances[detections_associees] = np.inf
        distances = sans_nan(distances)
        idx_best = int(np.argmin(distances))
        min_dist = distances[idx_best]

        if min_dist < seuil_mahalanobis:
            # CORRECTION
            mesure_Y = matrice_centres[idx_best]
            innovation = mesure_Y - C @ piste["X_k"]
            K = np.linalg.solve(S.T, (piste["P_k"] @ C.T).T).T     # P C' / S
            piste["X_k"] = piste["X_k"] + K @ innovation
            piste["P_k"] = (np.eye(9) - K @ C) @ piste["P_k"]

            piste["historique"] = ajouter_ligne(piste["historique"], piste["X_k"][:3])
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], mesure_Y)
            piste["lost_frames"] = 0
            piste["age"] += 1
            next_display_id = confirmer(piste, next_display_id, age_confirmation)
            detections_associees[idx_best] = True
        else:
            # OCCLUSION
            piste["lost_frames"] += 1
            piste["historique"] = ajouter_ligne(piste["historique"], piste["X_k"][:3])
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], [np.nan] * 3)

    # 3. NETTOYAGE
    tracks = [p for p in tracks if p["lost_frames"] <= trames_perdues_max]

    # 4. CRÉATION D'UNE NOUVELLE PISTE
    for d in range(nb_objets):
        if not detections_associees[d]:
            centre = matrice_centres[d]
            tracks.append({
                "id": next_id,
                "X_k": np.concatenate([centre, np.zeros(6)]),
                "P_k": np.eye(9) * P_initiale,
                "historique": centre.reshape(1, 3).copy(),
                "historique_observe": centre.reshape(1, 3).copy(),
                "lost_frames": 0,
                "age": 1,
                "confirmed": False,
                "display_id": 0,
            })
            next_id += 1

    return tracks, detections_associees, next_display_id, next_id
