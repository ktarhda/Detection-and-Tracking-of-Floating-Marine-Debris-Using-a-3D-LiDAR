"""Archive des vitesses (structure archive_vitesses de sauvegarde.m) et fichiers .mat.

Une entrée par objet (display_id), avec les mêmes champs que MATLAB. Chaque champ
(sauf display_id) est un vecteur qui grandit d'une valeur à chaque trame où la
piste est confirmée. Les fichiers .mat écrits ici s'ouvrent aussi dans MATLAB
(variable archive_vitesses_filtree, tableau de structures 1 x K), et les fichiers
.mat produits par MATLAB se relisent avec charger_mat.
"""

import numpy as np
import scipy.io

CHAMPS_ARCHIVE = (
    "display_id", "vitesses_moy", "vitesses_min", "vitesses_max",
    "x", "y", "z",
    "vx_moy", "vx_min", "vx_max",
    "vy_moy", "vy_min", "vy_max",
    "vz_moy", "vz_min", "vz_max",
    "ax_moy", "ax_min", "ax_max",
    "ay_moy", "ay_min", "ay_max",
    "az_moy", "az_min", "az_max",
    "p11", "p22", "p33", "p44", "p55", "p66", "p77", "p88", "p99",
    "frames", "est_observe",
)
CHAMPS_SERIES = CHAMPS_ARCHIVE[1:]
NOM_VARIABLE_MAT = "archive_vitesses_filtree"


class Archive:
    """archive_vitesses : liste d'entrées, dans l'ordre d'apparition des objets."""

    def __init__(self):
        self.entrees = []
        self._index = {}          # display_id -> position dans la liste

    def __len__(self):
        return len(self.entrees)

    def __iter__(self):
        return iter(self.entrees)

    def __getitem__(self, k):
        return self.entrees[k]

    def ajouter(self, display_id, etat, frame_idx, bool_observe):
        """Ajoute les valeurs d'une trame pour l'objet display_id (bloc « Variables pour l'archive »)."""
        P_diag = etat["P_diag"]
        pos = etat["pos_actuelle"]
        valeurs = {
            "vitesses_moy": etat["vitesses_moy"], "vitesses_min": etat["vitesses_min"],
            "vitesses_max": etat["vitesses_max"],
            "x": pos[0], "y": pos[1], "z": pos[2],
        }
        for axe in ("vx", "vy", "vz", "ax", "ay", "az"):
            for suffixe in ("moy", "min", "max"):
                valeurs[f"{axe}_{suffixe}"] = etat[f"{axe}_{suffixe}"]
        for k in range(9):
            valeurs[f"p{k + 1}{k + 1}"] = P_diag[k]
        valeurs["frames"] = frame_idx
        valeurs["est_observe"] = bool(bool_observe)

        position = self._index.get(display_id)
        if position is None:
            entree = {"display_id": display_id}
            for champ in CHAMPS_SERIES:
                entree[champ] = [valeurs[champ]]
            self._index[display_id] = len(self.entrees)
            self.entrees.append(entree)
        else:
            entree = self.entrees[position]
            for champ in CHAMPS_SERIES:
                entree[champ].append(valeurs[champ])

    def en_tableaux(self):
        """Copie de l'archive où chaque série est un tableau numpy."""
        return [_entree_en_tableaux(e) for e in self.entrees]


def _entree_en_tableaux(entree):
    sortie = {"display_id": int(entree["display_id"])}
    for champ in CHAMPS_SERIES:
        type_numpy = bool if champ == "est_observe" else float
        sortie[champ] = np.asarray(entree[champ], dtype=type_numpy).reshape(-1)
    return sortie


def filtrer_archive(archive, duree_min=10, ratio_min=0.9):
    """Section 7 de sauvegarde.m : on garde les objets vus assez longtemps et assez souvent.

    archive : liste d'entrées (tableaux numpy). Renvoie (archive_filtree, lignes_bilan)
    où lignes_bilan contient (display_id, duree_vie, ratio_observation, garde).
    """
    archive_filtree, bilan = [], []
    for entree in archive:
        duree_vie = len(entree["frames"])
        ratio_observation = float(np.sum(entree["est_observe"])) / duree_vie
        garde = duree_vie > duree_min and ratio_observation > ratio_min
        bilan.append((entree["display_id"], duree_vie, ratio_observation, garde))
        if garde:
            archive_filtree.append(entree)
    return archive_filtree, bilan


# ---------------------------------------------------------------------- fichiers .mat
def _tableau_structures(entrees):
    """Liste de dicts -> tableau de structures MATLAB 1 x K pour scipy.io.savemat."""
    type_struct = [(champ, object) for champ in CHAMPS_ARCHIVE]
    tableau = np.empty((1, len(entrees)), dtype=type_struct)
    for k, entree in enumerate(entrees):
        tableau[0, k]["display_id"] = float(entree["display_id"])
        for champ in CHAMPS_SERIES:
            valeurs = np.asarray(entree[champ])
            if champ == "est_observe":
                valeurs = valeurs.astype(bool)
            else:
                valeurs = valeurs.astype(float)
            tableau[0, k][champ] = valeurs.reshape(-1, 1)        # vecteur colonne, comme MATLAB
    return tableau


def sauvegarder_mat(chemin, archive_filtree, historique_cardinalite=None, n_reel=None):
    """save(nom_fichier, 'archive_vitesses_filtree', 'historique_cardinalite', 'n_reel')."""
    contenu = {NOM_VARIABLE_MAT: _tableau_structures(archive_filtree)}
    if historique_cardinalite is not None:
        contenu["historique_cardinalite"] = np.asarray(historique_cardinalite, dtype=float).reshape(-1, 1)
    if n_reel is not None:
        contenu["n_reel"] = float(n_reel)
    scipy.io.savemat(str(chemin), contenu, do_compression=True)


def _vers_dict(structure):
    """mat_struct (scipy) -> dict de tableaux numpy 1-D."""
    entree = {}
    for champ in structure._fieldnames:
        valeur = getattr(structure, champ)
        if champ == "display_id":
            entree[champ] = int(np.asarray(valeur).reshape(-1)[0])
        elif champ == "est_observe":
            entree[champ] = np.atleast_1d(np.asarray(valeur)).astype(bool).reshape(-1)
        else:
            entree[champ] = np.atleast_1d(np.asarray(valeur, dtype=float)).reshape(-1)
    return entree


def charger_mat(chemin):
    """load(nom_fichier) : renvoie un dict comme la structure de MATLAB.

    Clés : 'archive_vitesses_filtree' (liste de dicts), et si elles sont présentes
    'historique_cardinalite' (tableau 1-D) et 'n_reel' (entier).
    Fonctionne avec les fichiers écrits par sauvegarder_mat et avec ceux de MATLAB (format v7).
    """
    donnees = scipy.io.loadmat(str(chemin), squeeze_me=True, struct_as_record=False)
    if NOM_VARIABLE_MAT not in donnees:
        raise KeyError(f"La variable {NOM_VARIABLE_MAT} est absente de {chemin}.")
    brut = donnees[NOM_VARIABLE_MAT]
    if isinstance(brut, np.ndarray):
        structures = list(brut.reshape(-1))
    else:
        structures = [brut]
    sortie = {NOM_VARIABLE_MAT: [_vers_dict(s) for s in structures]}
    if "historique_cardinalite" in donnees:
        sortie["historique_cardinalite"] = np.atleast_1d(
            np.asarray(donnees["historique_cardinalite"], dtype=float)).reshape(-1)
    if "n_reel" in donnees:
        sortie["n_reel"] = int(np.asarray(donnees["n_reel"]).reshape(-1)[0])
    return sortie
