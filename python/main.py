"""main.m en Python : détection et suivi avec le filtre de Kalman classique ou ensembliste (UBIKF).

Comme main.m : affichage 3D en direct, bilan et figure de l'erreur de cardinalité,
puis graphiques de vitesse et d'accélération lissés (moyenne glissante sur 9 trames).
Ce script n'enregistre pas de fichier .mat (c'est le rôle de sauvegarde.py).

Exemples :
    python main.py
    python main.py --filtre ENSEMBLISTE
    python main.py --pcap 3objets.pcap --json 3objets.json --sans-visu
"""

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
