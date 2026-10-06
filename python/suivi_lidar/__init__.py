"""Détection et suivi de déchets flottants avec un LiDAR Ouster et une IMU.

Version Python du projet MATLAB (Lidar Toolbox, Image Processing Toolbox, INTLAB).

Modules :
    parametres              réglages du projet (valeurs des scripts MATLAB)
    lecture_ouster          ousterFileReader, readFrame, readIMU
    correction_inclinaison  correct_inclinaison.m
    roi                     région d'intérêt (main.m)
    segmentation            segmentation_rho_marin.m
    intervalles             arithmétique d'intervalles (INTLAB)
    filtres                 kalman_classique.m, kalman_ensembliste.m,
                            filtre_particule.m, filtre_particule_boite.m
    extraction              valeurs affichées et archivées (sauvegarde.m)
    archive                 archive_vitesses et fichiers .mat
    pipeline                boucle principale commune à main.m, Detect_particule.m et sauvegarde.m
    figures                 graphiques des scripts et de plot_figure.m
    visualisation           affichage 3D en direct (pcplayer)
    ligne_de_commande       options communes des scripts
    synthetique             enregistrement de test (scène simulée)
"""

from .parametres import (BOX_PARTICULAIRE, CLASSIQUE, ENSEMBLISTE, FILTRES, PARTICULAIRE, Parametres,
                         nom_filtre)

__version__ = "1.0.0"

__all__ = [
    "Parametres",
    "nom_filtre",
    "FILTRES",
    "CLASSIQUE",
    "ENSEMBLISTE",
    "PARTICULAIRE",
    "BOX_PARTICULAIRE",
]
