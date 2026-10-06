"""sauvegarde.m en Python : les quatre filtres, avec enregistrement des résultats.

Comme sauvegarde.m : affichage 3D en direct, bilan des pistes, fichier
donnees_<FILTRE>.mat (archive_vitesses_filtree, historique_cardinalite, n_reel)
et graphiques de vitesse et d'accélération. Les quatre fichiers .mat servent
ensuite à plot_figure.py.

Exemples :
    python sauvegarde.py --filtre CLASSIQUE
    python sauvegarde.py --filtre tous --sans-visu      (les quatre filtres à la suite)
"""

from suivi_lidar.ligne_de_commande import lancer_script

# =========================================================================
# CONFIGURATION DU FILTRE (utilisée quand --filtre n'est pas donné)
# =========================================================================
# CHOIX_FILTRE = "CLASSIQUE"
# CHOIX_FILTRE = "ENSEMBLISTE"
# CHOIX_FILTRE = "PARTICULAIRE"
CHOIX_FILTRE = "BOX_PARTICULAIRE"


def main(argv=None):
    return lancer_script("sauvegarde", CHOIX_FILTRE, argv,
                         description="Détection et suivi avec les quatre filtres et enregistrement des "
                                     "fichiers .mat (sauvegarde.m).")


if __name__ == "__main__":
    main()
