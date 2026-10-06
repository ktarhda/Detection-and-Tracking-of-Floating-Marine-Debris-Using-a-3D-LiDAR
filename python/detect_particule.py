

from suivi_lidar.ligne_de_commande import lancer_script

# =========================================================================
# SÉLECTION DU FILTRE (utilisée quand --filtre n'est pas donné)
# 'PARTICULAIRE' (nuage de points) ou 'BOX_PARTICULAIRE' (boîtes d'intervalles)
# =========================================================================
CHOIX_FILTRE = "BOX_PARTICULAIRE"
#CHOIX_FILTRE = "PARTICULAIRE"


def main(argv=None):
    return lancer_script("detect_particule", CHOIX_FILTRE, argv,
                         description="Détection et suivi : filtre particulaire ou Box Particle Filter "
                                     "(Detect_particule.m).")


if __name__ == "__main__":
    main()
