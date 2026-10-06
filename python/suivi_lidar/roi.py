"""Extraction de la région d'intérêt (partie ROI de main.m / sauvegarde.m)."""

import numpy as np


def colonnes_roi(n_colonnes, demi_angle=90.0, theta_avant=180.0):
    """Indices (à partir de 0) des colonnes face à l'eau.

    MATLAB : theta_brut = (0:N_cols-1) * (360/N_cols);
             ROI = abs(theta_brut - theta_avant) <= demi_angle;
    """
    theta = np.arange(n_colonnes) * (360.0 / n_colonnes)
    return np.flatnonzero(np.abs(theta - theta_avant) <= demi_angle)


def appliquer_roi(xyz_aligne, rho_brut, col_roi, portee_min=-0.4, portee_max=-15.0, limite_gauche=-7.0):
    """Garde les colonnes de la ROI et met à zéro les distances hors du plan d'eau.

    xyz_aligne : (lignes, colonnes, 3) nuage corrigé, NaN quand il n'y a pas de retour
    rho_brut   : (lignes, colonnes) distances en mètres, 0 quand il n'y a pas de retour
    Renvoie (rho_sans_sol, xyz_roi), comme rho_sans_sol et ptCloud_sans_sol dans MATLAB.
    """
    xyz_roi = xyz_aligne[:, col_roi, :]
    X, Y = xyz_roi[..., 0], xyz_roi[..., 1]
    with np.errstate(invalid="ignore"):
        zone_lac = (X >= portee_max) & (X <= portee_min) & (Y >= limite_gauche)
    rho_sans_sol = np.array(rho_brut[:, col_roi], dtype=float, copy=True)
    rho_sans_sol[~zone_lac] = 0.0
    return rho_sans_sol, xyz_roi
