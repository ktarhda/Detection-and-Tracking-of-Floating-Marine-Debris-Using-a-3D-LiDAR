"""Graphiques de main.m, Detect_particule.m, sauvegarde.m et plot_figure.m, avec matplotlib.

Chaque figure porte le même nom de fenêtre, le même titre, les mêmes couleurs,
les mêmes épaisseurs de trait et le même lissage (movmean) que dans MATLAB.
Elle est aussi enregistrée en PNG dans le dossier demandé.
"""

import re
import unicodedata
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np

from .parametres import BOX_PARTICULAIRE, CLASSIQUE, ENSEMBLISTE, PARTICULAIRE

# Couleurs MATLAB
ROUGE, VERT, BLEU, NOIR = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)
ORANGE = (1.0, 0.5, 0.0)
ORDRE_COULEURS_MATLAB = np.array([
    [0.0000, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250],
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330],
    [0.6350, 0.0780, 0.1840],
])
LARGEUR_DEFAUT = 0.5      # LineWidth par défaut de MATLAB

# Couleurs, épaisseurs et noms de plot_figure.m
STYLE_COMPARAISON = {
    CLASSIQUE: dict(couleur=ROUGE, largeur=1.5, nom="Filtre Classique", court="Classique"),
    ENSEMBLISTE: dict(couleur=VERT, largeur=2.0, nom="Filtre Ensembliste", court="Ensembliste"),
    PARTICULAIRE: dict(couleur=BLEU, largeur=1.5, nom="Filtre Particulaire", court="Particulaire"),
    BOX_PARTICULAIRE: dict(couleur=ORANGE, largeur=2.0, nom="Box Particulaire", court="Box Particulaire"),
}
ORDRE_FILTRES = (CLASSIQUE, ENSEMBLISTE, PARTICULAIRE, BOX_PARTICULAIRE)

# Les trois scripts MATLAB qui tracent les graphiques de vitesse
STYLE_MAIN = "main"
STYLE_DETECT_PARTICULE = "detect_particule"
STYLE_SAUVEGARDE = "sauvegarde"


def _plt():
    import matplotlib.pyplot as plt
    return plt


def movmean(x, k):
    """movmean(x, k) de MATLAB : moyenne glissante centrée, fenêtre réduite aux bords.

    k impair : (k-1)/2 valeurs de chaque côté ; k pair : k/2 avant et k/2 - 1 après.
    Un NaN dans la fenêtre donne NaN (option par défaut 'includenan').
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    n = x.size
    if n == 0 or k <= 1:
        return x.copy()
    avant = k // 2 if k % 2 == 0 else (k - 1) // 2
    apres = k // 2 - 1 if k % 2 == 0 else (k - 1) // 2
    sortie = np.empty(n)
    for i in range(n):
        sortie[i] = np.mean(x[max(0, i - avant):min(n, i + apres + 1)])
    return sortie


def hsv_matlab(n):
    """hsv(n) de MATLAB."""
    h = np.arange(n) / max(n, 1)
    return mcolors.hsv_to_rgb(np.column_stack([h, np.ones(n), np.ones(n)]))


def lines_matlab(n):
    """lines(n) de MATLAB."""
    return ORDRE_COULEURS_MATLAB[np.arange(n) % len(ORDRE_COULEURS_MATLAB)]


def _nom_fichier(nom):
    texte = unicodedata.normalize("NFKD", nom).encode("ascii", "ignore").decode("ascii")
    texte = re.sub(r"[^A-Za-z0-9]+", "_", texte).strip("_")
    return texte + ".png"


class GestionnaireFigures:
    """Crée les figures, les enregistre en PNG et les affiche à la fin si demandé."""

    def __init__(self, dossier=None, afficher=True):
        self.dossier = Path(dossier) if dossier is not None else None
        self.afficher = afficher
        self.fichiers = []
        if self.dossier is not None:
            self.dossier.mkdir(parents=True, exist_ok=True)

    def nouvelle(self, nom, **kwargs):
        plt = _plt()
        fig = plt.figure(**kwargs)
        fig.patch.set_facecolor("w")
        try:
            fig.canvas.manager.set_window_title(nom)
        except AttributeError:
            pass
        fig._nom_matlab = nom
        return fig

    def terminer(self, fig):
        if self.dossier is not None:
            chemin = self.dossier / _nom_fichier(fig._nom_matlab)
            fig.savefig(chemin, dpi=120, bbox_inches="tight")
            self.fichiers.append(chemin)
        if not self.afficher:
            _plt().close(fig)

    def montrer(self):
        if self.afficher:
            _plt().show()


def _axes_texte(ax, xlabel=None, ylabel=None, titre=None, taille_titre=12):
    if xlabel:
        ax.set_xlabel(xlabel, fontweight="bold")
    if ylabel:
        ax.set_ylabel(ylabel, fontweight="bold")
    if titre:
        ax.set_title(titre, fontsize=taille_titre)
    ax.grid(True)


# ====================================================================== graphiques de fin de script
# Différences entre les trois scripts MATLAB (version finale)
_TEXTES = {
    STYLE_SAUVEGARDE: dict(
        nom_g1="Avant Suppression des anomalies - {f}", titre_g1="Évolution de la vitesse globale (Avec Fantômes) - {f}",
        nom_g2="Après Suppression des anomalies - {f}", nom_g3="Décomposition Vitesse - ID {id} ({f})",
        titre_g3="Composantes 3D de la vitesse (Vx, Vy, Vz) - OBJET ID {id}", ylabel_acc="Accélération (m/s$^2$)",
        message_vide="Aucun objet n'a passé le filtre de ratio > 0.9. (Graphiques propres non générés)"),
    STYLE_MAIN: dict(
        nom_g1="Avant Supression d'anomalie (Tous les IDs)", titre_g1="Évolution de la vitesse globale (Avec Fantômes)",
        nom_g2="Après Supression d'anomalie - {f}", nom_g3="Décomposition Vitesse - ID {id}",
        titre_g3="Vitesses Estimées (Vx, Vy, Vz) - OBJET ID {id}", ylabel_acc="Accélération (m/s²)",
        message_vide=None),
    STYLE_DETECT_PARTICULE: dict(
        nom_g1="Avant Supression des anomalies - {f}", titre_g1="Évolution de la vitesse globale (Avec Fantômes) - {f}",
        nom_g2="Après Supression des anomalies - {f}", nom_g3="Décomposition Vitesse - ID {id}",
        titre_g3="Composantes 3D de la vitesse (Vx, Vy, Vz) - ID {id}", ylabel_acc="Accélération (m/s²)",
        message_vide="Aucun objet n'a passé le filtre de ratio > 0.5. (Graphiques propres non générés)"),
}


def figures_filtre(archive, archive_filtree, choix_filtre, gestionnaire, style=STYLE_SAUVEGARDE,
                   fenetre_lissage=9):
    """Graphiques 1 à 4 de fin de script (vitesse globale avant et après filtrage, Vx Vy Vz, ax ay az).

    style : 'sauvegarde' (sans lissage, couloirs pour les trois filtres à intervalles),
            'main' (lissage movmean, couloir seulement pour le filtre ensembliste),
            'detect_particule' (lissage movmean, couloirs toujours, borne basse >= 0).
    """
    if not archive:
        return
    t = _TEXTES[style]
    f_nom = choix_filtre
    if style == STYLE_SAUVEGARDE:
        avec_couloir = choix_filtre in (ENSEMBLISTE, PARTICULAIRE, BOX_PARTICULAIRE)
    elif style == STYLE_MAIN:
        avec_couloir = choix_filtre == ENSEMBLISTE
    else:
        avec_couloir = True

    def lisser(v):
        return np.asarray(v, dtype=float) if style == STYLE_SAUVEGARDE else movmean(v, fenetre_lissage)

    # GRAPHIQUE 1 : LA VITESSE GLOBALE (AVANT Suppression - Avec Fantômes)
    fig = gestionnaire.nouvelle(t["nom_g1"].format(f=f_nom))
    ax = fig.add_subplot(111)
    couleurs = hsv_matlab(len(archive))
    for k, entree in enumerate(archive):
        ax.plot(entree["frames"], lisser(entree["vitesses_moy"]), "-o", markersize=6, markerfacecolor="none",
                color=couleurs[k], linewidth=LARGEUR_DEFAUT, label=f"ID {entree['display_id']}")
    _axes_texte(ax, "Temps (Frames)", "Vitesse Absolue (m/s)", t["titre_g1"].format(f=f_nom))
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, ncol=2)
    gestionnaire.terminer(fig)

    if not archive_filtree:
        if t["message_vide"]:
            print(t["message_vide"])
        return

    # GRAPHIQUE 2 : LA VITESSE GLOBALE (APRÈS Suppression)
    fig = gestionnaire.nouvelle(t["nom_g2"].format(f=f_nom))
    ax = fig.add_subplot(111)
    couleurs = lines_matlab(len(archive_filtree))
    for k, entree in enumerate(archive_filtree):
        f = entree["frames"]
        v_moy = entree["vitesses_moy"]
        v_min = entree["vitesses_min"]
        if style == STYLE_DETECT_PARTICULE:
            v_min = np.maximum(0.0, v_min)
        moyenne_globale = float(np.mean(v_moy))
        etiquette = f"ID {entree['display_id']} (Moy: {moyenne_globale:.2f} m/s)"
        if avec_couloir:
            ax.fill_between(f, v_min, entree["vitesses_max"], color=couleurs[k], alpha=0.2, linewidth=0)
            ax.plot(f, lisser(v_moy), "-", linewidth=1.5, color=couleurs[k], label=etiquette)
        else:
            ax.plot(f, lisser(v_moy), "-", linewidth=LARGEUR_DEFAUT, color=couleurs[k], label=etiquette)
        ax.plot([f.min(), f.max()], [moyenne_globale, moyenne_globale], "--", linewidth=2, color=couleurs[k])
    _axes_texte(ax, "Temps (Frames)", "Vitesse Absolue (m/s)", f"Évolution de la vitesse globale ({f_nom})")
    ax.legend(loc="best", fontsize=10)
    gestionnaire.terminer(fig)

    # GRAPHIQUE 3 : DÉCOMPOSITION (Vx, Vy, Vz)
    for entree in archive_filtree:
        id_obj = entree["display_id"]
        f = entree["frames"]
        fig = gestionnaire.nouvelle(t["nom_g3"].format(id=id_obj, f=f_nom))
        ax = fig.add_subplot(111)
        if avec_couloir:
            for axe, couleur in (("vx", ROUGE), ("vy", VERT), ("vz", BLEU)):
                ax.fill_between(f, entree[f"{axe}_min"], entree[f"{axe}_max"], color=couleur, alpha=0.15,
                                linewidth=0)
        ax.plot(f, lisser(entree["vx_moy"]), "-", color=ROUGE, linewidth=1.5, label="Vx (Axe d'approche X)")
        ax.plot(f, lisser(entree["vy_moy"]), "-", color=VERT, linewidth=1.5, label="Vy (Axe latéral Y)")
        ax.plot(f, lisser(entree["vz_moy"]), "-", color=BLEU, linewidth=1.5, label="Vz (Axe vertical Z)")
        ax.plot([f.min(), f.max()], [0, 0], "--", color=NOIR, linewidth=1)
        _axes_texte(ax, "Temps (Frames)", "Vitesse (m/s)", t["titre_g3"].format(id=id_obj))
        ax.legend(loc="best", fontsize=10)
        gestionnaire.terminer(fig)

    # GRAPHIQUE 4 : DÉCOMPOSITION ACCÉLÉRATION (ax, ay, az)
    for entree in archive_filtree:
        id_obj = entree["display_id"]
        f = entree["frames"]
        fig = gestionnaire.nouvelle(f"Décomposition Accélération - ID {id_obj}")
        ax = fig.add_subplot(111)
        if avec_couloir:
            for axe, couleur in (("ax", ROUGE), ("ay", VERT), ("az", BLEU)):
                ax.fill_between(f, entree[f"{axe}_min"], entree[f"{axe}_max"], color=couleur, alpha=0.15,
                                linewidth=0)
        ax.plot(f, lisser(entree["ax_moy"]), "-", color=ROUGE, linewidth=1.5, label="ax (Axe X)")
        ax.plot(f, lisser(entree["ay_moy"]), "-", color=VERT, linewidth=1.5, label="ay (Axe Y)")
        ax.plot(f, lisser(entree["az_moy"]), "-", color=BLEU, linewidth=1.5, label="az (Axe Z)")
        ax.plot([f.min(), f.max()], [0, 0], "--", color=NOIR, linewidth=1)
        _axes_texte(ax, "Temps (Frames)", t["ylabel_acc"],
                    f"Accélérations Estimées (ax, ay, az) - OBJET ID {id_obj}")
        ax.legend(loc="best", fontsize=10)
        gestionnaire.terminer(fig)


def figure_cardinalite_filtre(historique_cardinalite, n_reel, choix_filtre, gestionnaire):
    """Figure « Erreur de Cardinalite - FILTRE » de main.m et Detect_particule.m."""
    card = np.asarray(historique_cardinalite, dtype=float)
    frames_card = np.arange(1, card.size + 1)
    erreur = np.abs(n_reel - card)
    fig = gestionnaire.nouvelle(f"Erreur de Cardinalite - {choix_filtre}")
    ax1, ax2 = fig.subplots(2, 1)
    ax1.plot(frames_card, n_reel * np.ones(card.size), "--", color=NOIR, linewidth=1.5,
             label="Nombre reel de cibles")
    ax1.plot(frames_card, card, "-", color=BLEU, linewidth=1.5, label="Nombre estime de cibles")
    ax1.grid(True)
    ax1.set_ylabel("Nombre de cibles", fontweight="bold")
    ax1.set_title(f"Estimation de la cardinalite - {choix_filtre}", fontsize=12)
    ax1.legend(loc="best")
    ax1.set_ylim(0, max(n_reel, card.max()) + 1)
    ax2.plot(frames_card, erreur, "-", color=ROUGE, linewidth=1.5)
    ax2.grid(True)
    ax2.set_xlabel("Temps (Frames)", fontweight="bold")
    ax2.set_ylabel("Erreur |n - m|", fontweight="bold")
    ax2.set_title("Erreur de cardinalite au cours du temps", fontsize=12)
    ax2.set_ylim(0, erreur.max() + 1)
    fig.tight_layout()
    gestionnaire.terminer(fig)


# ====================================================================== plot_figure.m
def _courbes(ax, donnees, champ, etiquettes=None, largeurs=None, fenetre=9):
    """Trace champ (lissé par movmean) pour les quatre filtres."""
    for filtre in ORDRE_FILTRES:
        _, entree, idx = donnees[filtre]
        style = STYLE_COMPARAISON[filtre]
        largeur = style["largeur"] if largeurs is None else largeurs
        etiquette = None if etiquettes is None else etiquettes[filtre]
        ax.plot(donnees["frames_communes"], movmean(entree[champ][idx], fenetre), "-", color=style["couleur"],
                linewidth=largeur, label=etiquette)


def figures_comparaison_objet(id_obj, donnees, gestionnaire, fenetre_lissage=9):
    """Les sept graphiques superposés de plot_figure.m pour un objet commun."""
    fen = fenetre_lissage
    courts = {flt: STYLE_COMPARAISON[flt]["court"] for flt in ORDRE_FILTRES}

    # GRAPHIQUE 1 : TRAJECTOIRE 3D SUPERPOSÉE
    fig = gestionnaire.nouvelle(f"Comparaison Trajectoire 3D - ID {id_obj}")
    ax = fig.add_subplot(111, projection="3d")
    for filtre in ORDRE_FILTRES:
        _, entree, idx = donnees[filtre]
        style = STYLE_COMPARAISON[filtre]
        ax.plot(movmean(entree["x"][idx], fen), movmean(entree["y"][idx], fen), movmean(entree["z"][idx], fen),
                "-", color=style["couleur"], linewidth=style["largeur"], label=style["nom"])
    _, entree_c, idx_c = donnees[CLASSIQUE]
    ax.plot([entree_c["x"][idx_c[0]]], [entree_c["y"][idx_c[0]]], [entree_c["z"][idx_c[0]]], "o",
            color=NOIR, markersize=8, markerfacecolor=NOIR, label="Départ")
    ax.view_init(elev=30, azim=-37.5)          # view(3)
    ax.set_xlabel("Position X (m)", fontweight="bold")
    ax.set_ylabel("Position Y (m)", fontweight="bold")
    ax.set_zlabel("Position Z (m)", fontweight="bold")
    ax.set_title(f"Comparaison des Trajectoires 3D - Objet ID {id_obj}", fontsize=12)
    ax.legend(loc="best")
    gestionnaire.terminer(fig)

    # GRAPHIQUE 2 : VITESSE GLOBALE
    fig = gestionnaire.nouvelle(f"Comparaison Vitesse Globale - ID {id_obj}")
    ax = fig.add_subplot(111)
    _courbes(ax, donnees, "vitesses_moy", {flt: STYLE_COMPARAISON[flt]["nom"] for flt in ORDRE_FILTRES},
             fenetre=fen)
    _axes_texte(ax, "Temps (Frames)", "Vitesse Absolue (m/s)",
                f"Comparaison de la Vitesse Globale Estimée - Objet ID {id_obj}")
    ax.legend(loc="best")
    gestionnaire.terminer(fig)

    # GRAPHIQUES À TROIS SOUS-FIGURES
    def trois_sous_figures(nom, champs, ylabels, titre, etiquettes, largeurs=None):
        fig = gestionnaire.nouvelle(nom, figsize=(7, 8))
        axes = fig.subplots(3, 1)
        for n, (ax, champ, ylabel) in enumerate(zip(axes, champs, ylabels)):
            _courbes(ax, donnees, champ, etiquettes if n == 0 else None, largeurs, fenetre=fen)
            ax.set_ylabel(ylabel, fontweight="bold")
            ax.grid(True)
        axes[0].set_title(titre)
        axes[0].legend(loc="best")
        axes[2].set_xlabel("Temps (Frames)", fontweight="bold")
        fig.tight_layout()
        gestionnaire.terminer(fig)

    # GRAPHIQUE 3 : COMPOSANTES DE LA VITESSE (Vx, Vy, Vz)
    trois_sous_figures(f"Vitesses Vx, Vy, Vz - ID {id_obj}", ("vx_moy", "vy_moy", "vz_moy"),
                       ("Vx (m/s)", "Vy (m/s)", "Vz (m/s)"), f"Décomposition Vitesse - ID {id_obj}",
                       courts, largeurs=1.5)

    # GRAPHIQUE 4 : INCERTITUDE DE POSITION (P11, P22, P33)
    trois_sous_figures(f"Incertitude Position - ID {id_obj}", ("p11", "p22", "p33"),
                       ("$P_{11}$ (X)", "$P_{22}$ (Y)", "$P_{33}$ (Z)"),
                       f"Comparaison des Covariances de Position - ID {id_obj}",
                       {CLASSIQUE: "Classique $P_{11}$", ENSEMBLISTE: "Ensembliste $P^+_{11}$",
                        PARTICULAIRE: "Particulaire Var(X)", BOX_PARTICULAIRE: "Box Particulaire (Largeur)"})

    # GRAPHIQUE 5 : INCERTITUDE DE VITESSE (P44, P55, P66)
    trois_sous_figures(f"Incertitude Vitesse - ID {id_obj}", ("p44", "p55", "p66"),
                       ("$P_{44}$ (Vx)", "$P_{55}$ (Vy)", "$P_{66}$ (Vz)"),
                       f"Comparaison des Covariances de Vitesse - ID {id_obj}",
                       {CLASSIQUE: "Classique $P_{44}$", ENSEMBLISTE: "Ensembliste $P^+_{44}$",
                        PARTICULAIRE: "Particulaire Var(Vx)", BOX_PARTICULAIRE: "Box Particulaire (Largeur)"})

    # GRAPHIQUE 6 : COMPOSANTES DE L'ACCÉLÉRATION (ax, ay, az)
    trois_sous_figures(f"Accélération ax, ay, az - ID {id_obj}", ("ax_moy", "ay_moy", "az_moy"),
                       ("ax (m/s$^2$)", "ay (m/s$^2$)", "az (m/s$^2$)"), f"Décomposition Accélération - ID {id_obj}",
                       courts, largeurs=1.5)

    # GRAPHIQUE 7 : INCERTITUDE D'ACCÉLÉRATION (P77, P88, P99)
    trois_sous_figures(f"Incertitude Accélération - ID {id_obj}", ("p77", "p88", "p99"),
                       ("$P_{77}$ (ax)", "$P_{88}$ (ay)", "$P_{99}$ (az)"),
                       f"Comparaison des Covariances d'Accélération - ID {id_obj}",
                       {CLASSIQUE: "Classique $P_{77}$", ENSEMBLISTE: "Ensembliste $P^+_{77}$",
                        PARTICULAIRE: "Particulaire Var(ax)", BOX_PARTICULAIRE: "Box Particulaire (Largeur)"})


def figure_cardinalite_comparaison(series, gestionnaire):
    """Figure « Comparaison Erreur de Cardinalite » de plot_figure.m.

    series : {filtre: (historique_cardinalite, n_reel)} pour les filtres qui ont ces données.
    Chaque filtre est tracé sur ses propres trames 1..n, comme dans MATLAB.
    """
    fig = gestionnaire.nouvelle("Comparaison Erreur de Cardinalite")
    ax1, ax2 = fig.subplots(2, 1)
    n_reel_ref = None
    for filtre in ORDRE_FILTRES:
        if filtre not in series:
            continue
        card, n_reel_ref = series[filtre]
        card = np.asarray(card, dtype=float)
        erreur = np.abs(n_reel_ref - card)
        fr = np.arange(1, card.size + 1)
        style = STYLE_COMPARAISON[filtre]
        ax1.plot(fr, card, "-", color=style["couleur"], linewidth=1.5, label=style["court"])
        ax2.plot(fr, erreur, "-", color=style["couleur"], linewidth=1.5, label=style["court"])
    if n_reel_ref is not None:
        limites = ax1.get_xlim()
        ax1.plot(limites, [n_reel_ref, n_reel_ref], "--", color=NOIR, linewidth=1.5, label="Nombre reel de cibles")
        ax1.set_xlim(limites)
    ax1.grid(True)
    ax1.set_ylabel("Nombre de cibles", fontweight="bold")
    ax1.set_title("Estimation de la cardinalite - Comparaison des 4 filtres", fontsize=12)
    ax1.legend(loc="best")
    ax2.grid(True)
    ax2.set_xlabel("Temps (Frames)", fontweight="bold")
    ax2.set_ylabel("Erreur |n - m|", fontweight="bold")
    ax2.set_title("Erreur de cardinalite au cours du temps", fontsize=12)
    ax2.legend(loc="best")
    fig.tight_layout()
    gestionnaire.terminer(fig)
