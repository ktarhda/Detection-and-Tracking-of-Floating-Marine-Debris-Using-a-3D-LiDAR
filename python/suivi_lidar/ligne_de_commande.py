"""Options de la ligne de commande, communes à main.py, detect_particule.py et sauvegarde.py."""

import argparse
import sys
from pathlib import Path

from .parametres import FILTRES, Parametres, nom_filtre


def preparer_console():
    """Évite une erreur d'encodage des symboles (✅, ❌, 💾) sur certaines consoles Windows.

    Un caractère que la console ne sait pas afficher est remplacé par « ? » au lieu
    d'arrêter le programme.
    """
    for flux in (sys.stdout, sys.stderr):
        try:
            flux.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


def creer_parser(description, choix_defaut, filtres_possibles, accepte_tous):
    defaut = Parametres()
    aide_filtre = ", ".join(filtres_possibles) + (", ou tous" if accepte_tous else "")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--filtre", default=choix_defaut, help=f"{aide_filtre} (défaut : %(default)s)")
    parser.add_argument("--pcap", default=defaut.fichier_pcap, help="fichier .pcap (défaut : %(default)s)")
    parser.add_argument("--json", default=defaut.fichier_json, help="fichier .json (défaut : %(default)s)")
    parser.add_argument("--premiere-trame", type=int, default=defaut.premiere_trame,
                        help="première trame traitée, comptée à partir de 1 (défaut : %(default)s)")
    parser.add_argument("--trames-ignorees-fin", type=int, default=defaut.trames_ignorees_fin,
                        help="nombre de trames ignorées à la fin du fichier (défaut : %(default)s)")
    parser.add_argument("--sortie", default="resultats",
                        help="dossier des fichiers .mat et des images (défaut : %(default)s)")
    parser.add_argument("--sans-visu", action="store_true", help="pas d'affichage 3D pendant le traitement")
    parser.add_argument("--sans-figures", action="store_true", help="pas de graphiques à la fin")
    parser.add_argument("--ne-pas-afficher", action="store_true",
                        help="enregistre les graphiques en PNG sans ouvrir de fenêtre")
    parser.add_argument("--graine", type=int, default=None,
                        help="graine du hasard des filtres particulaires (pour des résultats reproductibles)")
    parser.add_argument("--silencieux", action="store_true", help="moins de messages pendant le traitement")
    return parser


def lancer_script(nom_profil, choix_defaut, argv=None, description=""):
    """Lit les options puis lance le traitement pour le ou les filtres demandés."""
    from .pipeline import PROFILS

    preparer_console()
    profil = PROFILS[nom_profil]
    accepte_tous = set(profil.filtres) == set(FILTRES)
    args = creer_parser(description, choix_defaut, profil.filtres, accepte_tous).parse_args(argv)

    if accepte_tous and args.filtre.strip().lower() in ("tous", "all"):
        filtres = FILTRES
    else:
        try:
            filtres = (nom_filtre(args.filtre),)
        except ValueError as erreur:
            sys.exit(str(erreur))
        if filtres[0] not in profil.filtres:
            sys.exit(f"Ce script lance seulement : {', '.join(profil.filtres)}.")

    for fichier in (args.pcap, args.json):
        if not Path(fichier).is_file():
            sys.exit(f"Fichier introuvable : {fichier}. Donnez son chemin avec --pcap et --json.")

    # Avec plusieurs filtres à la suite, les graphiques sont seulement enregistrés
    # (sinon les fenêtres du premier filtre bloqueraient le suivant).
    afficher_figures = not args.ne_pas_afficher and len(filtres) == 1
    if not afficher_figures:
        import matplotlib
        matplotlib.use("Agg")

    from .pipeline import executer

    resultats = []
    for choix in filtres:
        params = Parametres(fichier_pcap=args.pcap, fichier_json=args.json,
                            premiere_trame=args.premiere_trame, trames_ignorees_fin=args.trames_ignorees_fin)
        resultats.append(executer(choix, params, profil=nom_profil, visualiser=not args.sans_visu,
                                  dossier_sortie=args.sortie, figures=not args.sans_figures,
                                  afficher_figures=afficher_figures, graine=args.graine,
                                  bavard=not args.silencieux))
        if len(filtres) > 1:
            print()
    return resultats
