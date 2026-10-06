"""Segmentation par distance radiale adaptative (segmentation_rho_marin.m).

Même algorithme que la fonction MATLAB :
1. balayage horizontal puis vertical de la matrice des distances ; un saut
   plus grand que le seuil marque une frontière ;
2. coupure de l'image des points valides le long des frontières ;
3. fermeture morphologique 3x3 (imclose) ;
4. étiquetage en 8-connexité (bwlabel), dans le même ordre que MATLAB ;
5. une boîte englobante par région d'au moins 5 points.

Règle du seuil, exactement comme dans le code MATLAB :
    - au début d'une ligne (ou après une coupure), le seuil vaut seuil_initial = 0,05 m
      et un saut est une frontière s'il dépasse ce seuil ;
    - ensuite, seuil = max(precision_ouster, moyenne des sauts de la surface en cours)
      et un saut est une frontière s'il dépasse 1,5 x seuil.
"""

import numpy as np
from scipy import ndimage

_CARRE_3x3 = np.ones((3, 3), dtype=bool)


def _balayage(rho, seuil_initial, precision, facteur):
    """Balayage ligne par ligne (gauche -> droite). Renvoie les frontières (lignes, colonnes-1)."""
    lignes, cols = rho.shape
    frontiere = np.zeros((lignes, cols - 1), dtype=bool)
    valides = (rho[:, 1:] > 0) & (rho[:, :-1] > 0)
    sauts = np.abs(rho[:, 1:] - rho[:, :-1])
    for i in range(lignes):
        js = np.flatnonzero(valides[i])
        if js.size == 0:
            continue
        seuil = seuil_initial
        somme, nombre = 0.0, 0
        sauts_i = sauts[i]
        for j in js:
            saut = sauts_i[j]
            limite = seuil if nombre == 0 else facteur * seuil
            if saut > limite:
                frontiere[i, j] = True
                seuil = seuil_initial
                somme, nombre = 0.0, 0
            else:
                somme += saut
                nombre += 1
                seuil = max(precision, somme / nombre)
    return frontiere


def segmentation_rho_marin(rho_sans_sol, xyz_roi, seuil_initial=0.05, precision_ouster=0.03,
                           facteur_seuil=1.5, nb_points_min=5):
    """Renvoie (candidats, boites).

    rho_sans_sol : (lignes, colonnes) distances en mètres, 0 = pas de mesure
    xyz_roi      : (lignes, colonnes, 3) coordonnées cartésiennes correspondantes
    candidats    : liste de dicts {'centre': array(3), 'nb_imp': int}
    boites       : tableau (n, 9) [x, y, z, largeur, profondeur, hauteur, 0, 0, 0]
    """
    rho = np.asarray(rho_sans_sol, dtype=float)

    # 1. Balayage horizontal (gauche -> droite) puis vertical (haut -> bas)
    frontiere_horiz = _balayage(rho, seuil_initial, precision_ouster, facteur_seuil)
    frontiere_vert = _balayage(rho.T, seuil_initial, precision_ouster, facteur_seuil).T

    # 2. Image des objets valides, coupée le long des frontières
    image = rho > 0
    image[:, :-1] &= ~frontiere_horiz
    image[:, 1:] &= ~frontiere_horiz
    image[:-1, :] &= ~frontiere_vert
    image[1:, :] &= ~frontiere_vert

    # 3. Fermeture morphologique 3x3 (comme imclose : bord à 0 pour la dilatation, à 1 pour l'érosion)
    dilatee = ndimage.binary_dilation(image, structure=_CARRE_3x3, border_value=0)
    image = ndimage.binary_erosion(dilatee, structure=_CARRE_3x3, border_value=1)

    # 4. Étiquetage 8-connexité. On étiquette la transposée pour retrouver l'ordre de
    #    bwlabel (parcours colonne par colonne), donc le même ordre des objets que MATLAB.
    labels_t, nb_regions = ndimage.label(image.T, structure=_CARRE_3x3)
    labels = labels_t.T

    # 5. Boîtes englobantes
    X, Y, Z = xyz_roi[..., 0], xyz_roi[..., 1], xyz_roi[..., 2]
    candidats, boites = [], []
    if nb_regions > 0:
        objets = ndimage.find_objects(labels)
        for k in range(1, nb_regions + 1):
            zone = objets[k - 1]
            masque = labels[zone] == k
            X_k, Y_k, Z_k = X[zone][masque], Y[zone][masque], Z[zone][masque]
            ok = ~np.isnan(X_k) & ~np.isnan(Y_k) & ~np.isnan(Z_k)
            X_k, Y_k, Z_k = X_k[ok], Y_k[ok], Z_k[ok]
            if X_k.size < nb_points_min:
                continue
            larg = X_k.max() - X_k.min()
            prof = Y_k.max() - Y_k.min()
            haut = Z_k.max() - Z_k.min()
            centre = np.array([X_k.min() + larg / 2, Y_k.min() + prof / 2, Z_k.min() + haut / 2])
            candidats.append({"centre": centre, "nb_imp": int(X_k.size)})
            boites.append([centre[0], centre[1], centre[2], larg, prof, haut, 0.0, 0.0, 0.0])

    boites = np.array(boites, dtype=float).reshape(-1, 9)
    return candidats, boites
