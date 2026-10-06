"""Extraction des valeurs affichées et archivées pour une piste confirmée.

Bloc « EXTRACTION SPÉCIFIQUE SELON LE FILTRE » de sauvegarde.m : pour chaque
filtre, on calcule la position, la vitesse et l'accélération (valeur centrale,
borne basse, borne haute), la diagonale de P et le texte affiché à côté de l'objet.
"""

import numpy as np

from .intervalles import Intervalle
from .parametres import BOX_PARTICULAIRE, CLASSIQUE, ENSEMBLISTE, PARTICULAIRE

AXES_VITESSE = ("vx", "vy", "vz")
AXES_ACCELERATION = ("ax", "ay", "az")


def _norme_intervalle(composantes):
    """sqrt(Vx^2 + Vy^2 + Vz^2) en arithmétique d'intervalles (comme INTLAB)."""
    total = composantes[0].carre() + composantes[1].carre() + composantes[2].carre()
    return total.racine()


def _valeurs_point(etat):
    """Cas où l'incertitude n'est pas un intervalle : min = max = moyenne."""
    v = {}
    for k, nom in enumerate(AXES_VITESSE):
        v[f"{nom}_moy"] = v[f"{nom}_min"] = v[f"{nom}_max"] = float(etat[3 + k])
    for k, nom in enumerate(AXES_ACCELERATION):
        v[f"{nom}_moy"] = v[f"{nom}_min"] = v[f"{nom}_max"] = float(etat[6 + k])
    return v


def _valeurs_intervalle(boite, moyennes=None):
    """Bornes inf/sup de chaque composante d'une boîte d'état (9 intervalles).

    moyennes : valeurs centrales à utiliser (sinon mid de chaque intervalle).
    """
    if moyennes is None:
        moyennes = boite.mid()
    v = {}
    for k, nom in enumerate(AXES_VITESSE + AXES_ACCELERATION):
        i = 3 + k
        v[f"{nom}_moy"] = float(moyennes[i])
        v[f"{nom}_min"] = float(boite.inf[i])
        v[f"{nom}_max"] = float(boite.sup[i])
    return v


def extraire_etat(piste, choix_filtre, bpf_enveloppe=False):
    """Renvoie un dict avec toutes les valeurs de l'archive et de l'affichage.

    Clés : pos_actuelle (3,), vitesses_moy, vitesses_min, vitesses_max,
    vx_moy ... az_max, P_diag (9,), texte_affichage,
    et selon le filtre : boite_position (UBIKF) ou boites (BPF), pour le dessin.

    bpf_enveloppe : pour le BPF, valeurs centrales prises au milieu de l'enveloppe
    des boîtes (Detect_particule.m) au lieu de la moyenne pondérée X_k_estime
    (sauvegarde.m).
    """
    display_id = piste["display_id"]
    sortie = {}

    if choix_filtre == CLASSIQUE:
        X_k = piste["X_k"]
        vitesse_moy = float(np.linalg.norm(X_k[3:6]))
        sortie.update(_valeurs_point(X_k))
        sortie.update(pos_actuelle=X_k[0:3].copy(), vitesses_moy=vitesse_moy,
                      vitesses_min=vitesse_moy, vitesses_max=vitesse_moy,
                      P_diag=np.diag(piste["P_k"]).copy(),
                      texte_affichage=f"ID {display_id} | {vitesse_moy:.2f} m/s")

    elif choix_filtre == ENSEMBLISTE:
        X = piste["X_k_intval"]
        V_norm_int = _norme_intervalle([X[3], X[4], X[5]])
        milieux = X.mid()
        vitesse_moy = float(np.sqrt(milieux[3] ** 2 + milieux[4] ** 2 + milieux[5] ** 2))
        v_min, v_max = float(V_norm_int.inf), float(V_norm_int.sup)
        sortie.update(_valeurs_intervalle(X))
        sortie.update(pos_actuelle=milieux[0:3].copy(), vitesses_moy=vitesse_moy,
                      vitesses_min=v_min, vitesses_max=v_max,
                      P_diag=np.diag(piste["P_k_plus"]).copy(),
                      texte_affichage=f"ID {display_id} | V:[{v_min:.2f}, {v_max:.2f}]",
                      boite_position=X[0:3].copy())

    elif choix_filtre == PARTICULAIRE:
        X_k = piste["X_k"]
        particules = piste["particules"]
        vitesse_moy = float(np.linalg.norm(X_k[3:6]))
        ecarts = particules.std(axis=1, ddof=1)               # std(particules(k,:)) de MATLAB
        for k, nom in enumerate(AXES_VITESSE + AXES_ACCELERATION):
            moy = float(X_k[3 + k])
            sortie[f"{nom}_moy"] = moy
            sortie[f"{nom}_min"] = moy - 2 * float(ecarts[3 + k])
            sortie[f"{nom}_max"] = moy + 2 * float(ecarts[3 + k])
        marge = 2 * float(np.linalg.norm(ecarts[3:6]))
        sortie.update(pos_actuelle=X_k[0:3].copy(), vitesses_moy=vitesse_moy,
                      vitesses_min=max(0.0, vitesse_moy - marge), vitesses_max=vitesse_moy + marge,
                      P_diag=particules.var(axis=1, ddof=1),     # var(particules, 0, 2)
                      texte_affichage=f"ID {display_id} | {vitesse_moy:.2f} m/s")

    elif choix_filtre == BOX_PARTICULAIRE:
        boites = piste["X_boxes"]
        # Enveloppe de toutes les boîtes (min des bornes inf, max des bornes sup ; NaN ignorés comme min/max de MATLAB)
        X_hull = Intervalle.infsup(np.fmin.reduce(boites.inf, axis=1), np.fmax.reduce(boites.sup, axis=1))
        V_norm_int = _norme_intervalle([X_hull[3], X_hull[4], X_hull[5]])
        v_min, v_max = float(V_norm_int.inf), float(V_norm_int.sup)
        if bpf_enveloppe:
            # Detect_particule.m : valeurs centrales = milieu de l'enveloppe
            milieux = X_hull.mid()
            vitesse_moy = float(np.sqrt(milieux[3] ** 2 + milieux[4] ** 2 + milieux[5] ** 2))
            sortie.update(_valeurs_intervalle(X_hull))
            position = milieux[0:3].copy()
        else:
            # sauvegarde.m : valeurs centrales = moyenne pondérée des boîtes (X_k_estime)
            X_estime = piste["X_k_estime"]
            vitesse_moy = float(np.linalg.norm(X_estime[3:6]))
            sortie.update(_valeurs_intervalle(X_hull, moyennes=X_estime))
            position = X_estime[0:3].copy()
        sortie.update(pos_actuelle=position, vitesses_moy=vitesse_moy,
                      vitesses_min=v_min, vitesses_max=v_max,
                      P_diag=X_hull.rad(),
                      texte_affichage=f"ID {display_id} | V:[{v_min:.2f}, {v_max:.2f}]",
                      boites=boites, poids=piste["poids"])
    else:
        raise ValueError(f"Filtre inconnu : {choix_filtre}")

    return sortie
