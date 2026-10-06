"""Filtre particulaire multi-objets (filtre_particule.m).

Champs d'une piste : id, X_k (9,), particules (9, N), poids (N,), historique,
historique_observe, lost_frames, age, confirmed, display_id.
"""

import numpy as np

from ..parametres import matrice_A, matrice_C
from .commun import ajouter_ligne, confirmer, distances_mahalanobis, nb_detections


def filtre_particule(tracks, matrice_centres, dt, N_particules, Q_cov, R_cov, next_display_id, next_id,
                     rng=None, seuil_mahalanobis=3.5, age_confirmation=5, trames_perdues_max=10,
                     ecarts_types_initiaux=None, ecart_max=None):
    rng = np.random.default_rng() if rng is None else rng
    nb_objets = nb_detections(matrice_centres)
    detections_associees = np.zeros(nb_objets, dtype=bool)
    if ecarts_types_initiaux is None:
        ecarts_types_initiaux = np.array([0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1])
    if ecart_max is None:
        ecart_max = np.array([3, 3, 3, 2, 2, 2, 2, 2, 2], dtype=float)   # dispersion max (pos, vit, acc)
    ecart_max = np.asarray(ecart_max, dtype=float)[:, None]

    # Matrices du modèle cinématique à 9 états
    A = matrice_A(dt)
    C = matrice_C()
    inv_R = np.linalg.inv(R_cov)
    L_Q = np.linalg.cholesky(Q_cov)        # chol(Q_cov)' de MATLAB (triangulaire inférieure)

    # 1. PRÉDICTION (déplacement de toutes les particules + bruit gaussien)
    for piste in tracks:
        bruit = L_Q @ rng.standard_normal((9, N_particules))
        piste["particules"] = A @ piste["particules"] + bruit
        piste["X_k"] = piste["particules"] @ piste["poids"]          # moyenne pondérée
        # Limitation de la dispersion autour de la moyenne
        ecarts = piste["particules"] - piste["X_k"][:, None]
        ecarts = np.maximum(np.minimum(ecarts, ecart_max), -ecart_max)
        piste["particules"] = piste["X_k"][:, None] + ecarts

    # 2. ASSOCIATION (MAHALANOBIS) ET CORRECTION
    for piste in tracks:
        if nb_objets == 0:
            break

        centre_predit = piste["X_k"][:3]
        P_nuage = np.cov(piste["particules"][:3], ddof=1)             # cov(particules(1:3,:)')
        S = P_nuage + R_cov
        distances = distances_mahalanobis(matrice_centres, centre_predit, np.linalg.inv(S),
                                          detections_associees)
        idx_best = int(np.argmin(distances))
        min_mahal = distances[idx_best]

        if min_mahal < seuil_mahalanobis:
            mesure_Y = matrice_centres[idx_best]

            # Mise à jour des poids par la vraisemblance gaussienne
            innov = mesure_Y[:, None] - C @ piste["particules"]
            dist2 = np.sum(innov * (inv_R @ innov), axis=0)
            vraisemblance = np.exp(-0.5 * dist2) + 1e-15
            poids = piste["poids"] * vraisemblance
            poids = poids / poids.sum()

            # Rééchantillonnage si N_eff < N/2 (randsample avec remise)
            N_eff = 1.0 / np.sum(poids ** 2)
            if N_eff < N_particules / 2:
                indices = rng.choice(N_particules, size=N_particules, replace=True, p=poids)
                piste["particules"] = piste["particules"][:, indices]
                poids = np.full(N_particules, 1.0 / N_particules)
            piste["poids"] = poids

            piste["X_k"] = piste["particules"] @ piste["poids"]
            piste["historique"] = ajouter_ligne(piste["historique"], piste["X_k"][:3])
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], mesure_Y)
            piste["lost_frames"] = 0
            piste["age"] += 1
            next_display_id = confirmer(piste, next_display_id, age_confirmation)
            detections_associees[idx_best] = True
        else:
            # COASTING (objet perdu)
            piste["lost_frames"] += 1
            piste["historique"] = ajouter_ligne(piste["historique"], piste["X_k"][:3])
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], [np.nan] * 3)

    # 3. NETTOYAGE DES FANTÔMES
    tracks = [p for p in tracks if p["lost_frames"] <= trames_perdues_max]

    # 4. CRÉATION D'UNE NOUVELLE PISTE
    L_init = np.linalg.cholesky(np.diag(np.asarray(ecarts_types_initiaux, dtype=float) ** 2))
    for d in range(nb_objets):
        if not detections_associees[d]:
            centre = matrice_centres[d]
            centre_initial = np.concatenate([centre, np.zeros(6)])
            particules = centre_initial[:, None] + L_init @ rng.standard_normal((9, N_particules))
            tracks.append({
                "id": next_id,
                "X_k": centre_initial,
                "particules": particules,
                "poids": np.full(N_particules, 1.0 / N_particules),
                "historique": centre.reshape(1, 3).copy(),
                "historique_observe": centre.reshape(1, 3).copy(),
                "lost_frames": 0,
                "age": 1,
                "confirmed": False,
                "display_id": 0,
            })
            next_id += 1

    return tracks, detections_associees, next_display_id, next_id
