

from suivi_lidar.ligne_de_commande import lancer_script

# =========================================================================
# CONFIGURATION DU FILTRE (utilisée quand --filtre n'est pas donné)
# =========================================================================
#CHOIX_FILTRE = "CLASSIQUE"
CHOIX_FILTRE = "ENSEMBLISTE"


def main(argv=None):
    return lancer_script("main", CHOIX_FILTRE, argv,
                         description="Détection et suivi : Kalman classique ou ensembliste (main.m).")


if __name__ == "__main__":
    main()
