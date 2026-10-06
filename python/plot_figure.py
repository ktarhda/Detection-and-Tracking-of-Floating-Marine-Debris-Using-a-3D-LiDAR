"""ÉTUDE COMPARATIVE : CLASSIQUE vs ENSEMBLISTE vs PARTICULAIRE vs BOX PARTICULAIRE (plot_figure.m).

Charge les quatre fichiers donnees_<FILTRE>.mat (écrits par sauvegarde.py, ou par
sauvegarde.m dans MATLAB), cherche les objets communs aux quatre filtres,
synchronise les trames communes et trace les sept graphiques superposés par objet
(courbes lissées par une moyenne glissante sur 9 trames), puis l'erreur de cardinalité.

Exemples :
    python plot_figure.py
    python plot_figure.py --dossier resultats --ne-pas-afficher
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from suivi_lidar.archive import NOM_VARIABLE_MAT, charger_mat
from suivi_lidar.ligne_de_commande import preparer_console
from suivi_lidar.parametres import BOX_PARTICULAIRE, CLASSIQUE, ENSEMBLISTE, PARTICULAIRE

fenetre_lissage = 9   # Taille de la fenêtre du filtre moyenne glissante (en frames)

ORDRE = (CLASSIQUE, ENSEMBLISTE, PARTICULAIRE, BOX_PARTICULAIRE)
NOMS_CHARGEMENT = {CLASSIQUE: "Classique", ENSEMBLISTE: "Ensembliste", PARTICULAIRE: "Particulaire",
                   BOX_PARTICULAIRE: "Box Particulaire"}


def synchroniser(frames_par_filtre):
    """Synchronisation temporelle stricte, exactement comme plot_figure.m.

    frames_par_filtre : liste de quatre vecteurs de trames (C, E, P, BP).
    Renvoie (frames_communes, [idx_C, idx_E, idx_P, idx_BP]).
    """
    f_C, f_E, f_P, f_BP = (np.asarray(f) for f in frames_par_filtre)
    frames_temp1, iC_t1, iE_t1 = np.intersect1d(f_C, f_E, return_indices=True)
    frames_temp2, iCE_t2, iP_t2 = np.intersect1d(frames_temp1, f_P, return_indices=True)
    frames_communes, iCEP_t3, idx_BP = np.intersect1d(frames_temp2, f_BP, return_indices=True)
    idx_C = iC_t1[iCE_t2[iCEP_t3]]
    idx_E = iE_t1[iCE_t2[iCEP_t3]]
    idx_P = iP_t2[iCEP_t3]
    return frames_communes, [idx_C, idx_E, idx_P, idx_BP]


def lire_arguments(argv=None):
    parser = argparse.ArgumentParser(description="Comparaison des quatre filtres (plot_figure.m).")
    parser.add_argument("--dossier", default="resultats", help="dossier des fichiers .mat (défaut : %(default)s)")
    parser.add_argument("--figures", default=None,
                        help="dossier des images PNG (défaut : <dossier>/figures_comparaison)")
    parser.add_argument("--ne-pas-afficher", action="store_true", help="enregistre les PNG sans ouvrir de fenêtre")
    return parser.parse_args(argv)


def main(argv=None):
    preparer_console()
    args = lire_arguments(argv)
    if args.ne_pas_afficher:
        import matplotlib
        matplotlib.use("Agg")
    from suivi_lidar.figures import GestionnaireFigures, figure_cardinalite_comparaison, figures_comparaison_objet
    from suivi_lidar.pipeline import bilan_cardinalite

    print("  ÉTUDE COMPARATIVE : CLASSIQUE vs ENSEMBLISTE vs PARTICULAIRE vs BOX PARTICULAIRE  ")
    dossier = Path(args.dossier)

    # =========================================================================
    # 1. CHARGEMENT DES 4 FICHIERS .MAT
    # =========================================================================
    data = {}
    for filtre in ORDRE:
        chemin = dossier / f"donnees_{filtre}.mat"
        if not chemin.is_file():
            sys.exit(f"Fichier {chemin} introuvable. Lancez sauvegarde.py avec le filtre {filtre} d'abord.")
        data[filtre] = charger_mat(chemin)
        print(f"✅ Fichier {NOMS_CHARGEMENT[filtre]} chargé.")
    archives = {f: data[f][NOM_VARIABLE_MAT] for f in ORDRE}

    # =========================================================================
    # 2. RECHERCHE DES OBJETS COMMUNS (aux 4 méthodes)
    # =========================================================================
    ids = {f: np.array([e["display_id"] for e in archives[f]]) for f in ORDRE}
    ids_communs = np.intersect1d(np.intersect1d(np.intersect1d(ids[CLASSIQUE], ids[ENSEMBLISTE]),
                                                ids[PARTICULAIRE]), ids[BOX_PARTICULAIRE])
    if ids_communs.size == 0:
        print("❌ Aucun ID d'objet commun trouvé entre les quatre exécutions.")
        return
    print(f"🎯 {ids_communs.size} objet(s) commun(s) trouvé(s) pour la comparaison !")

    gestionnaire = GestionnaireFigures(Path(args.figures) if args.figures else dossier / "figures_comparaison",
                                       afficher=not args.ne_pas_afficher)

    # =========================================================================
    # 3. GÉNÉRATION DES GRAPHIQUES SUPERPOSÉS
    # =========================================================================
    for id_obj in ids_communs:
        objets = {f: archives[f][int(np.flatnonzero(ids[f] == id_obj)[0])] for f in ORDRE}
        # --- SYNCHRONISATION TEMPORELLE STRICTE (4 filtres) ---
        frames_communes, indices = synchroniser([objets[f]["frames"] for f in ORDRE])
        if frames_communes.size == 0:
            continue
        donnees = {"frames_communes": frames_communes}
        for filtre, idx in zip(ORDRE, indices):
            donnees[filtre] = (objets[filtre]["frames"], objets[filtre], idx)
        figures_comparaison_objet(int(id_obj), donnees, gestionnaire, fenetre_lissage)

    print("✅ Graphiques comparatifs des 4 filtres générés avec succès !")

    # =========================================================================
    # 4. MÉTRIQUE D'ÉVALUATION MOT : ERREUR DE CARDINALITÉ
    # =========================================================================
    print("\n================== ERREUR DE CARDINALITE ==================")
    print(f"{'Filtre':<18} | {'Moy est':>8} | {'Err moy':>8} | {'Exact (%)':>9} | {'Max':>5}")
    print("-----------------------------------------------------------")
    series = {}
    for filtre in ORDRE:
        if "historique_cardinalite" not in data[filtre]:
            print(f"{NOMS_CHARGEMENT[filtre]:<18} | (donnees absentes : relancer sauvegarde.py)")
            continue
        card = data[filtre]["historique_cardinalite"]
        n_reel = data[filtre]["n_reel"]
        b = bilan_cardinalite(card, n_reel)
        print(f"{NOMS_CHARGEMENT[filtre]:<18} | {b['moyenne']:8.2f} | {b['erreur_moyenne']:8.2f} | "
              f"{b['trames_exactes']:9.1f} | {b['erreur_max']:5d}")
        series[filtre] = (card, n_reel)
    print("===========================================================\n")
    figure_cardinalite_comparaison(series, gestionnaire)

    if gestionnaire.fichiers:
        print(f"-> {len(gestionnaire.fichiers)} images enregistrées dans {gestionnaire.dossier}")
    gestionnaire.montrer()


if __name__ == "__main__":
    main()
