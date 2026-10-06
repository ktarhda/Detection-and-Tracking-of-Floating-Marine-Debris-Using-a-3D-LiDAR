"""Petites fonctions communes aux quatre filtres."""

import numpy as np


def nb_detections(matrice_centres):
    """size(matrice_centres, 1) : 0 quand il n'y a aucune détection."""
    if matrice_centres is None:
        return 0
    m = np.asarray(matrice_centres)
    return 0 if m.size == 0 else m.shape[0]


def ajouter_ligne(historique, ligne):
    """historique = [historique; ligne] (ajout d'une ligne de 3 valeurs)."""
    return np.vstack([historique, np.asarray(ligne, dtype=float).reshape(1, 3)])


def confirmer(piste, next_display_id, age_confirmation):
    """Logique de confirmation commune aux quatre fonctions MATLAB.

    age vaut 1 à la création et augmente de 1 à chaque association. Quand
    age >= age_confirmation (5 ou 10 selon le filtre), la piste est confirmée
    et reçoit son identifiant affiché (display_id), une seule fois.
    """
    if piste["age"] >= age_confirmation:
        if not piste["confirmed"]:
            piste["display_id"] = next_display_id
            next_display_id += 1
        piste["confirmed"] = True
    return next_display_id


def distances_mahalanobis(matrice_centres, centre_predit, inv_S, detections_associees):
    """Distance de Mahalanobis entre chaque détection et la position prédite.

    Les détections déjà prises par une autre piste reçoivent Inf.
    """
    innov = matrice_centres - centre_predit
    d = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", innov, inv_S, innov), 0.0))
    d[detections_associees] = np.inf
    return sans_nan(d)


def sans_nan(distances):
    """min() de MATLAB ignore les NaN : on les remplace par Inf avant np.argmin."""
    distances = np.asarray(distances, dtype=float)
    distances[np.isnan(distances)] = np.inf
    return distances
