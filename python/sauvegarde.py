

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
