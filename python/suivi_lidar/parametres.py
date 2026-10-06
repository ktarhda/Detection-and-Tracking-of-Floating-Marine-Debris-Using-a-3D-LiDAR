"""Tous les paramètres du projet, au même endroit.

Les valeurs par défaut sont celles de la version finale du code MATLAB
(sauvegarde.m et les quatre fonctions de filtre). Les seules valeurs propres
à un script sont changées par ce script (par exemple erreur_lidar_max = 0.05
pour le filtre ensembliste dans main.m). Le commentaire de chaque champ
indique d'où vient la valeur.
"""

from dataclasses import dataclass, field

import numpy as np

# Noms des filtres, identiques à la variable CHOIX_FILTRE de MATLAB
CLASSIQUE = "CLASSIQUE"
ENSEMBLISTE = "ENSEMBLISTE"
PARTICULAIRE = "PARTICULAIRE"
BOX_PARTICULAIRE = "BOX_PARTICULAIRE"
FILTRES = (CLASSIQUE, ENSEMBLISTE, PARTICULAIRE, BOX_PARTICULAIRE)

# Noms courts acceptés en ligne de commande
ALIAS_FILTRES = {
    "classique": CLASSIQUE, "kalman": CLASSIQUE,
    "ensembliste": ENSEMBLISTE, "ubikf": ENSEMBLISTE,
    "particulaire": PARTICULAIRE, "pf": PARTICULAIRE,
    "box_particulaire": BOX_PARTICULAIRE, "box": BOX_PARTICULAIRE, "bpf": BOX_PARTICULAIRE,
}


def nom_filtre(texte):
    """Convertit un nom saisi par l'utilisateur en nom officiel de filtre."""
    cle = texte.strip()
    if cle.upper() in FILTRES:
        return cle.upper()
    if cle.lower() in ALIAS_FILTRES:
        return ALIAS_FILTRES[cle.lower()]
    raise ValueError(f"Filtre inconnu : {texte}. Choix possibles : {', '.join(FILTRES)}")


@dataclass
class Parametres:
    # ------------------------------------------------------------------ acquisition
    fichier_pcap: str = "3objets.pcap"
    fichier_json: str = "3objets.json"
    premiere_trame: int = 720            # for i = 720 : (total_frames - 470), indices à partir de 1
    trames_ignorees_fin: int = 470
    frequence_hz: int = 10
    dt: float = 0.1

    # ------------------------------------------------------------------ région d'intérêt
    demi_angle: float = 90.0
    theta_avant: float = 180.0
    portee_min: float = -0.4
    portee_max: float = -15.0
    limite_gauche: float = -7.0

    # ------------------------------------------------------------------ segmentation (segmentation_rho_marin.m)
    seuil_initial: float = 0.05
    precision_ouster: float = 0.03
    facteur_seuil: float = 1.5
    nb_points_min: int = 5

    # ------------------------------------------------------------------ association et vie des pistes
    seuil_mahalanobis: float = 3.5
    # Confirmation : « if tracks(t).age >= ... » dans chaque fonction de filtre
    age_confirmation: dict = field(default_factory=lambda: {
        CLASSIQUE: 5,            # kalman_classique.m
        ENSEMBLISTE: 10,         # kalman_ensembliste.m
        PARTICULAIRE: 5,         # filtre_particule.m
        BOX_PARTICULAIRE: 10,    # filtre_particule_boite.m
    })
    # Suppression : « pistes_valides = [tracks.lost_frames] <= 10 » dans les quatre fonctions
    trames_perdues_max: dict = field(default_factory=lambda: {
        CLASSIQUE: 10, ENSEMBLISTE: 10, PARTICULAIRE: 10, BOX_PARTICULAIRE: 10,
    })

    # ------------------------------------------------------------------ Kalman classique
    Q_classique: np.ndarray = field(default_factory=lambda: np.eye(9) * 0.01)
    R_classique: np.ndarray = field(default_factory=lambda: np.eye(3) * 0.02)
    P_initiale_classique: float = 10.0

    # ------------------------------------------------------------------ Kalman ensembliste / UBIKF
    erreur_lidar_max_ensembliste: float = 0.2    # sauvegarde.m (0.05 dans main.m)
    Q_ensembliste: np.ndarray = field(default_factory=lambda: np.eye(9) * 0.01)
    R_ensembliste: np.ndarray = field(default_factory=lambda: np.eye(3) * 0.02)
    q_bruit_jerk: float = 0.01                   # Q_opt de kalman_ensembliste.m
    facteur_sigma_contraction: float = 1.5

    # ------------------------------------------------------------------ bornes du modèle (V_modele, UBIKF et BPF)
    rayons_V_modele: np.ndarray = field(default_factory=lambda: np.array(
        [0.02, 0.02, 0.02, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5]))

    # ------------------------------------------------------------------ filtre particulaire
    N_particules: int = 1000
    Q_particulaire: np.ndarray = field(default_factory=lambda: np.diag(
        np.array([0.02, 0.02, 0.02, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05]) ** 2))
    R_particulaire: np.ndarray = field(default_factory=lambda: np.diag(np.array([0.03, 0.03, 0.03]) ** 2))
    ecarts_types_initiaux_particules: np.ndarray = field(default_factory=lambda: np.array(
        [0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1]))
    ecart_max_particules: np.ndarray = field(default_factory=lambda: np.array(
        [3, 3, 3, 2, 2, 2, 2, 2, 2], dtype=float))      # dispersion max autorisée (filtre_particule.m)

    # ------------------------------------------------------------------ Box Particle Filter
    N_boites: int = 10
    erreur_lidar_max_boites: float = 0.03
    R_boites: np.ndarray = field(default_factory=lambda: np.diag(np.array([0.03, 0.03, 0.03]) ** 2))
    rayon_max_boites: np.ndarray = field(default_factory=lambda: np.array(
        [2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5], dtype=float))
    rayon_max_association: float = 0.3

    # ------------------------------------------------------------------ incertitude des nouvelles pistes (UBIKF et BPF)
    rayons_initiaux_vitesse_acceleration: np.ndarray = field(default_factory=lambda: np.array(
        [0.5, 0.5, 0.5, 0.1, 0.1, 0.1]))

    # ------------------------------------------------------------------ archive, cardinalité et graphiques
    duree_min_archive: int = 10           # duree_vie > 10
    ratio_observation_min: float = 0.9    # ratio_observation > 0.9
    duree_validation: int = 50            # message « ID validé » si au moins 50 trames (sauvegarde.m)
    distance_max_couleur: float = 40.0
    n_reel: int = 3                       # nombre réel de cibles déployées dans le canal
    fenetre_lissage: int = 9              # moyenne glissante des graphiques (main.m, Detect_particule.m, plot_figure.m)

    # ------------------------------------------------------------------ matrices du modèle cinématique
    @property
    def A(self):
        return matrice_A(self.dt)

    @property
    def C(self):
        return matrice_C()


def matrice_A(dt):
    """Modèle à accélération constante, 9 états : position, vitesse, accélération."""
    I3, Z3 = np.eye(3), np.zeros((3, 3))
    return np.block([
        [I3, I3 * dt, I3 * (0.5 * dt ** 2)],
        [Z3, I3, I3 * dt],
        [Z3, Z3, I3],
    ])


def matrice_C():
    """Le LiDAR ne mesure que la position : C = [I3 0 0]."""
    return np.hstack([np.eye(3), np.zeros((3, 6))])


def matrice_Q_jerk(dt, q):
    """Q_opt de kalman_ensembliste.m : bruit blanc sur le jerk intégré sur dt."""
    I3 = np.eye(3)
    return np.block([
        [dt ** 5 / 20 * I3, dt ** 4 / 8 * I3, dt ** 3 / 6 * I3],
        [dt ** 4 / 8 * I3, dt ** 3 / 3 * I3, dt ** 2 / 2 * I3],
        [dt ** 3 / 6 * I3, dt ** 2 / 2 * I3, dt * I3],
    ]) * q
