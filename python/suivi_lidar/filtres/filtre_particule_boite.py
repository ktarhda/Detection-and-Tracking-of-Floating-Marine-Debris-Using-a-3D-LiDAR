"""Box Particle Filter multi-objets (filtre_particule_boite.m).

Chaque piste porte N boîtes (tableau d'intervalles 9 x N) et leurs poids.
Champs : id, X_boxes, X_prec_estime (9,), poids (N,), X_k_estime (9,),
historique, historique_observe, lost_frames, age, confirmed, display_id.
"""

import numpy as np

from ..intervalles import Intervalle, intersect, produit
from ..parametres import matrice_A, matrice_C
from .commun import ajouter_ligne, confirmer, distances_mahalanobis, nb_detections


def filtre_particule_boite(tracks, matrice_centres, dt, N_boites, V_modele, erreur_lidar_max, R_cov,
                           next_display_id, next_id, rng=None, seuil_mahalanobis=3.5, age_confirmation=10,
                           trames_perdues_max=10, rayon_max=None, rayon_gating_max=0.3,
                           rayons_initiaux_va=None):
    rng = np.random.default_rng() if rng is None else rng
    nb_objets = nb_detections(matrice_centres)
    detections_associees = np.zeros(nb_objets, dtype=bool)
    if rayon_max is None:
        rayon_max = np.array([2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5], dtype=float)
    if rayons_initiaux_va is None:
        rayons_initiaux_va = np.array([0.5, 0.5, 0.5, 0.1, 0.1, 0.1])

    A = matrice_A(dt)
    C = matrice_C()
    V = Intervalle(V_modele.inf[:, None], V_modele.sup[:, None])      # même bruit pour chaque boîte

    # 1. PRÉDICTION DES BOÎTES
    for piste in tracks:
        piste["X_prec_estime"] = piste["X_k_estime"].copy()
        boites = produit(A, piste["X_boxes"]) + V
        # Saturation des boîtes (évite l'explosion par wrapping effect)
        rayons = np.minimum(boites.rad(), rayon_max[:, None])
        piste["X_boxes"] = Intervalle.midrad(boites.mid(), rayons)
        centres = piste["X_boxes"][0:3].mid()
        piste["X_k_estime"][:3] = centres @ piste["poids"]

    # 2. ASSOCIATION ET CORRECTION
    for piste in tracks:
        if nb_objets == 0:
            break

        centre_predit = piste["X_k_estime"][:3]
        # Incertitude = dispersion des centres + largeur des boîtes (plafonnée pour la porte)
        centres = piste["X_boxes"][0:3].mid()
        rayons = np.minimum(piste["X_boxes"][0:3].rad(), rayon_gating_max)
        P_nuage = np.cov(centres, ddof=1) + np.diag(rayons.mean(axis=1) ** 2)
        S_k = P_nuage + R_cov
        distances = distances_mahalanobis(matrice_centres, centre_predit, np.linalg.inv(S_k),
                                          detections_associees)
        idx_best = int(np.argmin(distances))
        min_dist = distances[idx_best]

        if min_dist < seuil_mahalanobis:
            mesure_Y = matrice_centres[idx_best]
            boite_lidar = Intervalle.midrad(mesure_Y, erreur_lidar_max)
            poids_temp = np.zeros(N_boites)

            for b in range(N_boites):
                boite = piste["X_boxes"][:, b]
                Z_box = produit(C, boite)
                R_box = intersect(Z_box, boite_lidar)
                if R_box.contient_nan():
                    poids_temp[b] = 0.0
                    continue
                volume_Z = np.prod(Z_box.rad() * 2 + 1e-6)
                volume_R = np.prod(R_box.rad() * 2 + 1e-6)
                poids_temp[b] = piste["poids"][b] * (volume_R / volume_Z)

                # 1. Contraction de la position
                piste["X_boxes"][0:3, b] = R_box
                # 2. Contraction de la vitesse (backward)
                V_backward = (R_box - piste["X_prec_estime"][:3]) / dt
                V_contractee = intersect(piste["X_boxes"][3:6, b], V_backward)
                piste["X_boxes"][3:6, b] = V_backward if V_contractee.contient_nan() else V_contractee

            somme_poids = poids_temp.sum()
            if somme_poids > 0:
                piste["poids"] = poids_temp / somme_poids
                piste["lost_frames"] = 0
                N_eff = 1.0 / np.sum(piste["poids"] ** 2)
                if N_eff < N_boites / 2:
                    piste["X_boxes"] = subdivision_resampling(piste["X_boxes"], piste["poids"], N_boites, rng)
                    piste["poids"] = np.full(N_boites, 1.0 / N_boites)
            else:
                # Sécurité anti-fantômes : la porte a déjà validé l'association (min_dist < 3.5).
                # Position = boîte de mesure, vitesse = contrainte backward (accélération inchangée).
                V_backward = (boite_lidar - piste["X_prec_estime"][:3]) / dt
                for b in range(N_boites):
                    piste["X_boxes"][0:3, b] = boite_lidar
                    piste["X_boxes"][3:6, b] = V_backward
                piste["poids"] = np.full(N_boites, 1.0 / N_boites)
                piste["lost_frames"] = 0

            piste["X_k_estime"] = piste["X_boxes"].mid() @ piste["poids"]
            piste["historique"] = ajouter_ligne(piste["historique"], piste["X_k_estime"][:3])
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], mesure_Y)
            piste["age"] += 1
            next_display_id = confirmer(piste, next_display_id, age_confirmation)
            detections_associees[idx_best] = True
        else:
            piste["lost_frames"] += 1
            piste["historique"] = ajouter_ligne(piste["historique"], piste["X_k_estime"][:3])
            piste["historique_observe"] = ajouter_ligne(piste["historique_observe"], [np.nan] * 3)

    # 3. NETTOYAGE
    tracks = [p for p in tracks if p["lost_frames"] <= trames_perdues_max]

    # 4. CRÉATION D'UNE NOUVELLE PISTE
    for d in range(nb_objets):
        if not detections_associees[d]:
            centre = matrice_centres[d]
            centre_initial = np.concatenate([centre, np.zeros(6)])
            incertitude = np.concatenate([np.full(3, erreur_lidar_max), rayons_initiaux_va])
            boite_mere = Intervalle.midrad(centre_initial, incertitude)
            tracks.append({
                "id": next_id,
                "X_boxes": diviser_boite(boite_mere, N_boites, 3),   # dimension 4 en MATLAB (vx)
                "X_prec_estime": centre_initial.copy(),
                "poids": np.full(N_boites, 1.0 / N_boites),
                "X_k_estime": centre_initial.copy(),
                "historique": centre.reshape(1, 3).copy(),
                "historique_observe": centre.reshape(1, 3).copy(),
                "lost_frames": 0,
                "age": 1,
                "confirmed": False,
                "display_id": 0,
            })
            next_id += 1

    return tracks, detections_associees, next_display_id, next_id


def subdivision_resampling(boxes, weights, N, rng):
    """Rééchantillonnage systématique, puis découpe des boîtes choisies."""
    cum_w = np.cumsum(weights)
    step = 1.0 / N
    u = rng.random() * step
    counts = np.zeros(N, dtype=int)
    idx = 0
    for _ in range(N):
        while u > cum_w[idx] and idx < N - 1:
            idx += 1
        counts[idx] += 1
        u += step

    nouvelles = Intervalle.zeros((boxes.shape[0], N))
    new_idx = 0
    for i in range(N):
        k = counts[i]
        if k > 0:
            dim_max = int(np.argmax(boxes[:, i].rad()))        # plus grand côté
            sous = diviser_boite(boxes[:, i], k, dim_max)
            nouvelles[:, new_idx:new_idx + k] = sous
            new_idx += k
    return nouvelles


def diviser_boite(boite, k, dim):
    """Découpe une boîte en k sous-boîtes égales le long de la dimension dim (index à partir de 0)."""
    inf = np.repeat(boite.inf[:, None], k, axis=1)
    sup = np.repeat(boite.sup[:, None], k, axis=1)
    if k > 1:
        bas, haut = boite.inf[dim], boite.sup[dim]
        pas = (haut - bas) / k
        i = np.arange(k)
        inf[dim] = bas + i * pas
        sup[dim] = bas + (i + 1) * pas
    return Intervalle(inf, sup)
