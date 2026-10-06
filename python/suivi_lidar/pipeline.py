"""Boucle principale de détection et de suivi, commune à main.m, Detect_particule.m et sauvegarde.m.

Pour chaque trame : lecture, correction de l'inclinaison, ROI, segmentation,
appel du filtre choisi, cardinalité, extraction des valeurs des pistes confirmées,
archive, affichage. À la fin : bilans, fichier .mat (sauvegarde.m) et graphiques.

Les trois scripts MATLAB font presque la même chose ; leurs petites différences
(messages, réglages, lissage des graphiques, moment où la cardinalité est comptée...)
sont décrites dans un « profil » (voir PROFILS) pour que chaque script Python
reproduise exactement son script MATLAB.
"""

import time
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .archive import Archive, filtrer_archive, sauvegarder_mat
from .correction_inclinaison import correct_inclinaison, moyenne_accelerometre
from .extraction import extraire_etat
from .figures import (STYLE_DETECT_PARTICULE, STYLE_MAIN, STYLE_SAUVEGARDE, GestionnaireFigures,
                      figure_cardinalite_filtre, figures_filtre)
from .filtres import filtre_particule, filtre_particule_boite, kalman_classique, kalman_ensembliste
from .intervalles import Intervalle
from .lecture_ouster import LecteurOuster, resolution_json
from .parametres import (BOX_PARTICULAIRE, CLASSIQUE, ENSEMBLISTE, FILTRES, PARTICULAIRE, Parametres,
                         nom_filtre)
from .roi import appliquer_roi, colonnes_roi
from .segmentation import segmentation_rho_marin


@dataclass(frozen=True)
class ProfilScript:
    """Ce qui distingue main.m, Detect_particule.m et sauvegarde.m (version finale)."""
    nom: str
    filtres: tuple                       # filtres que le script sait lancer
    titre: str                           # premier message affiché
    filtres_intervalles: tuple           # filtres pour lesquels le script démarre INTLAB
    affiche_frequence: bool              # « -> Fréquence LiDAR détectée ... » (main.m)
    affiche_detections_par_trame: bool   # « Objets détectés sur la frame ... » (sauvegarde.m)
    cardinalite_avant_filtre: bool       # main.m compte les pistes confirmées avant l'appel du filtre
    bpf_enveloppe: bool                  # Detect_particule.m : valeurs du BPF au milieu de l'enveloppe
    trace_observe: bool                  # trajectoire observée en jaune (main.m, sauvegarde.m)
    bilan_pistes_finales: bool           # « ID validé / rejeté » sur les pistes encore actives (sauvegarde.m)
    rapport_cardinalite: bool            # bilan + figure de cardinalité (main.m, Detect_particule.m)
    bilan_archive: str                   # messages du filtrage de l'archive : 'main', 'detect_particule' ou ''
    enregistre_mat: bool                 # fichier donnees_<FILTRE>.mat (sauvegarde.m)
    message_aucun_objet: bool            # « Aucun objet n'a été détecté durant la session. »
    style_figures: str
    erreur_lidar_max_ensembliste: float = None    # valeur propre au script (main.m : 0.05)


PROFILS = {
    "main": ProfilScript(
        nom="main", filtres=(CLASSIQUE, ENSEMBLISTE),
        titre="  PERCEPTION MARINE : DÉTECTION (MODÈLE 9 ÉTATS)  ",
        filtres_intervalles=(ENSEMBLISTE,), affiche_frequence=True, affiche_detections_par_trame=False,
        cardinalite_avant_filtre=True, bpf_enveloppe=False, trace_observe=True, bilan_pistes_finales=False,
        rapport_cardinalite=True, bilan_archive="main", enregistre_mat=False, message_aucun_objet=False,
        style_figures=STYLE_MAIN, erreur_lidar_max_ensembliste=0.05),
    "detect_particule": ProfilScript(
        nom="detect_particule", filtres=(PARTICULAIRE, BOX_PARTICULAIRE),
        titre="  PERCEPTION MARINE : FILTRES PARTICULAIRES (Classique & Boîtes) ",
        filtres_intervalles=(BOX_PARTICULAIRE,), affiche_frequence=False, affiche_detections_par_trame=False,
        cardinalite_avant_filtre=False, bpf_enveloppe=True, trace_observe=False, bilan_pistes_finales=False,
        rapport_cardinalite=True, bilan_archive="detect_particule", enregistre_mat=False,
        message_aucun_objet=True, style_figures=STYLE_DETECT_PARTICULE),
    "sauvegarde": ProfilScript(
        nom="sauvegarde", filtres=FILTRES,
        titre="  PERCEPTION MARINE : DÉTECTION ET SUIVI  ",
        filtres_intervalles=(ENSEMBLISTE, BOX_PARTICULAIRE), affiche_frequence=False,
        affiche_detections_par_trame=True, cardinalite_avant_filtre=False, bpf_enveloppe=False,
        trace_observe=True, bilan_pistes_finales=True, rapport_cardinalite=False, bilan_archive="",
        enregistre_mat=True, message_aucun_objet=True, style_figures=STYLE_SAUVEGARDE),
}


class Suivi:
    """Variables du tracker (tracks, next_id, next_display_id) et appel du filtre choisi.

    Équivalent du bloc « VARIABLES SPÉCIFIQUES SELON LE FILTRE » et du bloc
    « APPEL DU TRACKER » des scripts MATLAB.
    """

    def __init__(self, choix_filtre, params=None, rng=None):
        self.choix_filtre = nom_filtre(choix_filtre)
        self.p = params if params is not None else Parametres()
        self.rng = rng if rng is not None else np.random.default_rng()
        self.dt = self.p.dt
        self.A = self.p.A
        self.C = self.p.C
        self.tracks = []
        self.next_id = 1
        self.next_display_id = 1
        if self.choix_filtre in (ENSEMBLISTE, BOX_PARTICULAIRE):
            # V_modele = midrad(zeros(9,1), [0.02; 0.02; 0.02; 0.1; 0.1; 0.1; 0.5; 0.5; 0.5])
            self.V_modele = Intervalle.midrad(np.zeros(9), self.p.rayons_V_modele)

    def mettre_a_jour(self, matrice_centres):
        """Un appel de la fonction du filtre ; renvoie detections_associees."""
        p, choix = self.p, self.choix_filtre
        matrice_centres = np.asarray(matrice_centres, dtype=float).reshape(-1, 3)
        reglages = dict(seuil_mahalanobis=p.seuil_mahalanobis, age_confirmation=p.age_confirmation[choix],
                        trames_perdues_max=p.trames_perdues_max[choix])
        if choix == CLASSIQUE:
            self.tracks, detections, self.next_display_id, self.next_id = kalman_classique(
                self.tracks, matrice_centres, self.dt, self.A, self.C, p.Q_classique, p.R_classique,
                self.next_display_id, self.next_id, P_initiale=p.P_initiale_classique, **reglages)
        elif choix == ENSEMBLISTE:
            self.tracks, detections, self.next_display_id, self.next_id = kalman_ensembliste(
                self.tracks, matrice_centres, self.dt, self.A, self.C, p.Q_ensembliste, p.R_ensembliste,
                self.V_modele, p.erreur_lidar_max_ensembliste, self.next_display_id, self.next_id,
                q_bruit=p.q_bruit_jerk, facteur_sigma=p.facteur_sigma_contraction,
                rayons_initiaux_va=p.rayons_initiaux_vitesse_acceleration, **reglages)
        elif choix == PARTICULAIRE:
            self.tracks, detections, self.next_display_id, self.next_id = filtre_particule(
                self.tracks, matrice_centres, self.dt, p.N_particules, p.Q_particulaire, p.R_particulaire,
                self.next_display_id, self.next_id, rng=self.rng,
                ecarts_types_initiaux=p.ecarts_types_initiaux_particules, ecart_max=p.ecart_max_particules,
                **reglages)
        else:
            self.tracks, detections, self.next_display_id, self.next_id = filtre_particule_boite(
                self.tracks, matrice_centres, self.dt, p.N_boites, self.V_modele, p.erreur_lidar_max_boites,
                p.R_boites, self.next_display_id, self.next_id, rng=self.rng, rayon_max=p.rayon_max_boites,
                rayon_gating_max=p.rayon_max_association,
                rayons_initiaux_va=p.rayons_initiaux_vitesse_acceleration, **reglages)
        return detections

    def nb_pistes_confirmees(self):
        """m_estime = sum([tracks.confirmed])."""
        return int(sum(1 for piste in self.tracks if piste["confirmed"]))

    def pistes_confirmees(self):
        return [piste for piste in self.tracks if piste["confirmed"]]


def bilan_cardinalite(historique_cardinalite, n_reel):
    """Erreur de cardinalité : moyenne estimée, erreur moyenne, % de trames exactes, erreur max."""
    card = np.asarray(historique_cardinalite, dtype=float)
    if card.size == 0:
        return {"moyenne": np.nan, "erreur_moyenne": np.nan, "trames_exactes": np.nan, "erreur_max": np.nan}
    erreur = np.abs(n_reel - card)
    return {
        "moyenne": float(card.mean()),
        "erreur_moyenne": float(erreur.mean()),
        "trames_exactes": float(100.0 * np.sum(erreur == 0) / erreur.size),
        "erreur_max": int(erreur.max()),
    }


def _afficher_rapport_cardinalite(historique_cardinalite, n_reel, choix):
    b = bilan_cardinalite(historique_cardinalite, n_reel)
    print(f"\n--- ERREUR DE CARDINALITE ({choix}) ---")
    print(f"Nombre reel de cibles          : {n_reel}")
    print(f"Nombre estime moyen            : {b['moyenne']:.2f}")
    print(f"Erreur de cardinalite moyenne  : {b['erreur_moyenne']:.2f}")
    print(f"Trames a cardinalite exacte    : {b['trames_exactes']:.1f} %")
    print(f"Erreur de cardinalite maximale : {b['erreur_max']}")
    print("--------------------------------------------------")


def executer(choix_filtre, params=None, profil="sauvegarde", visualiser=True, dossier_sortie="resultats",
             figures=True, afficher_figures=True, graine=None, bavard=True):
    """Lance la détection et le suivi sur tout l'enregistrement, pour un filtre.

    choix_filtre : 'CLASSIQUE', 'ENSEMBLISTE', 'PARTICULAIRE' ou 'BOX_PARTICULAIRE'
    params       : Parametres (fichiers, seuils, réglages des filtres)
    profil       : 'main', 'detect_particule' ou 'sauvegarde' (script MATLAB reproduit)
    visualiser   : affichage 3D en direct (Open3D)
    figures      : graphiques de fin de traitement (enregistrés en PNG, affichés si afficher_figures)
    graine       : graine du hasard (filtres particulaires) ; None = tirage différent à chaque fois, comme MATLAB
    Renvoie un dict : archive, archive_filtree, historique_cardinalite, n_reel, tracks, fichier_mat,
    figures, nb_trames.
    """
    profil = PROFILS[profil] if isinstance(profil, str) else profil
    p = params if params is not None else Parametres()
    choix = nom_filtre(choix_filtre)
    if choix not in profil.filtres:
        raise ValueError(f"Le script {profil.nom} ne lance que les filtres : {', '.join(profil.filtres)}.")
    if profil.erreur_lidar_max_ensembliste is not None:
        p = replace(p, erreur_lidar_max_ensembliste=profil.erreur_lidar_max_ensembliste)
    dossier_sortie = Path(dossier_sortie)
    dossier_sortie.mkdir(parents=True, exist_ok=True)

    print(profil.titre)
    print(f"-> Filtre : {choix}")
    if choix in profil.filtres_intervalles:
        print("-> Arithmétique d'intervalles : module intervalles.py (à la place d'INTLAB)")

    # LECTURE JSON + RÉSOLUTION, ROI
    N_cols, N_rows = resolution_json(p.fichier_json)
    col_ROI = colonnes_roi(N_cols, p.demi_angle, p.theta_avant)

    # LECTEUR + IMU
    lecteur = LecteurOuster(p.fichier_pcap, p.fichier_json)
    imu_data = lecteur.lire_imu()
    total_frames = lecteur.nombre_trames
    if imu_data["accelerometre"].shape[0] == 0:
        print("⚠️  Aucune mesure IMU dans le fichier : le nuage n'est pas corrigé (roll = pitch = 0).")
        accel_moyen = np.array([0.0, 0.0, 9.80665])
    else:
        accel_moyen = moyenne_accelerometre(imu_data)       # mean(AccelerometerReadings), calculée une fois
    if bavard:
        print(f"-> Fichier : {p.fichier_pcap} | {total_frames} trames | {N_rows} x {N_cols} | "
              f"{imu_data['accelerometre'].shape[0]} mesures IMU")
    if profil.affiche_frequence:
        print(f"-> Fréquence LiDAR détectée : {p.frequence_hz} Hz | dt = {p.dt:.3f} seconde")

    # VISUALISATEUR
    viewer = None
    if visualiser:
        try:
            from .visualisation import Visualiseur
            viewer = Visualiseur(choix, trace_observe=profil.trace_observe)
        except (ImportError, RuntimeError) as erreur:
            print(f"⚠️  Affichage 3D désactivé : {erreur}")
            viewer = None

    suivi = Suivi(choix, p, np.random.default_rng(graine))
    archive = Archive()
    historique_cardinalite = []

    # BOUCLE PRINCIPALE : for i = 720 : (total_frames - 470)
    derniere = total_frames - p.trames_ignorees_fin
    if derniere < p.premiere_trame:
        print(f"⚠️  Aucune trame à traiter : le fichier a {total_frames} trames, il en faut au moins "
              f"{p.premiere_trame + p.trames_ignorees_fin} avec premiere_trame = {p.premiere_trame} "
              f"et trames_ignorees_fin = {p.trames_ignorees_fin} (options --premiere-trame et "
              f"--trames-ignorees-fin).")
    frame_idx = 0
    debut = time.perf_counter()
    try:
        for i, xyz_brut, rho_brut in lecteur.trames(p.premiere_trame, derniere):
            if viewer is not None and not viewer.est_ouvert():
                break
            frame_idx += 1

            # 1. ACQUISITION + CALIBRATION
            xyz_aligne = correct_inclinaison(xyz_brut, accel_moyen)

            # 2. APPLICATION DU ROI
            rho_sans_sol, xyz_roi = appliquer_roi(xyz_aligne, rho_brut, col_ROI, p.portee_min, p.portee_max,
                                                  p.limite_gauche)

            # 3. SEGMENTATION DIRECTE
            candidats, boites = segmentation_rho_marin(rho_sans_sol, xyz_roi, p.seuil_initial,
                                                       p.precision_ouster, p.facteur_seuil, p.nb_points_min)
            nb_objets = len(candidats)
            if profil.affiche_detections_par_trame and bavard:
                print(f"Objets détectés sur la frame {frame_idx} : {nb_objets}")

            # 4. APPEL DU TRACKER (et métrique de cardinalité)
            matrice_centres = np.array([c["centre"] for c in candidats], dtype=float).reshape(-1, 3)
            if profil.cardinalite_avant_filtre:
                historique_cardinalite.append(suivi.nb_pistes_confirmees())
            suivi.mettre_a_jour(matrice_centres)
            if not profil.cardinalite_avant_filtre:
                historique_cardinalite.append(suivi.nb_pistes_confirmees())

            # 5. EXTRACTION ET ARCHIVE (pistes confirmées)
            pistes_affichees = []
            for piste in suivi.tracks:
                if piste["confirmed"]:
                    etat = extraire_etat(piste, choix, bpf_enveloppe=profil.bpf_enveloppe)
                    archive.ajouter(piste["display_id"], etat, frame_idx, piste["lost_frames"] == 0)
                    pistes_affichees.append((piste, etat))

            # 6. AFFICHAGE 3D
            if viewer is not None:
                from .visualisation import couleurs_nuage, points_affichables
                couleurs = couleurs_nuage(rho_brut, xyz_aligne, boites, p.distance_max_couleur)
                xyz_visu, couleurs_visu = points_affichables(xyz_aligne, couleurs)
                viewer.afficher(xyz_visu, couleurs_visu, boites, pistes_affichees)
    finally:
        lecteur.fermer()
        if viewer is not None:
            viewer.fermer()

    duree = time.perf_counter() - debut
    print("Fin de l'acquisition marine.")
    if bavard and frame_idx > 0:
        print(f"-> {frame_idx} trames traitées en {duree:.1f} s ({frame_idx / max(duree, 1e-9):.1f} trames/s)")

    historique_cardinalite = np.asarray(historique_cardinalite, dtype=float)
    resultats = {"filtre": choix, "profil": profil.nom, "archive": [], "archive_filtree": [],
                 "historique_cardinalite": historique_cardinalite, "n_reel": p.n_reel, "tracks": suivi.tracks,
                 "fichier_mat": None, "figures": [], "nb_trames": frame_idx}

    # BILAN DES PISTES ENCORE ACTIVES (sauvegarde.m)
    if profil.bilan_pistes_finales:
        print("--------------------------------------------------")
        for piste in suivi.tracks:
            duree_piste = piste["historique"].shape[0]
            if piste["display_id"] > 0:
                if duree_piste >= p.duree_validation:
                    print(f"✅ ID {piste['display_id']} validé : Durée = {duree_piste} frames")
                else:
                    print(f"❌ ID {piste['display_id']} rejeté (fantôme/vague) : Durée = {duree_piste} frames")
        print("--------------------------------------------------")

    gestionnaire = None
    if figures:
        gestionnaire = GestionnaireFigures(dossier_sortie / f"figures_{profil.nom}_{choix}", afficher=afficher_figures)

    # MÉTRIQUE D'ÉVALUATION MOT : ERREUR DE CARDINALITÉ (main.m, Detect_particule.m)
    if profil.rapport_cardinalite and historique_cardinalite.size > 0:
        _afficher_rapport_cardinalite(historique_cardinalite, p.n_reel, choix)
        if gestionnaire is not None:
            figure_cardinalite_filtre(historique_cardinalite, p.n_reel, choix, gestionnaire)

    # FILTRAGE FINAL DE L'ARCHIVE ET SAUVEGARDE DU FICHIER .MAT
    archive_tableaux = archive.en_tableaux()
    resultats["archive"] = archive_tableaux
    if archive_tableaux:
        archive_filtree, bilan = filtrer_archive(archive_tableaux, p.duree_min_archive, p.ratio_observation_min)
        resultats["archive_filtree"] = archive_filtree
        for display_id, duree_vie, ratio, garde in bilan:
            if profil.bilan_archive == "main":
                if garde:
                    print(f"✅ ID {display_id} validé : Durée = {duree_vie} frames")
                else:
                    print(f"❌ ID {display_id} rejeté (fantôme/vague) : Durée = {duree_vie} frames")
            elif profil.bilan_archive == "detect_particule":
                if garde:
                    print(f"✅ ID {display_id} validé : Durée = {duree_vie} frames | Ratio = {ratio:.2f}")
                else:
                    print(f"❌ ID {display_id} rejeté (Fantôme) : Durée = {duree_vie} frames | Ratio = {ratio:.2f}")
        if profil.enregistre_mat and archive_filtree:
            nom_fichier = dossier_sortie / f"donnees_{choix}.mat"
            sauvegarder_mat(nom_fichier, archive_filtree, historique_cardinalite, p.n_reel)
            resultats["fichier_mat"] = nom_fichier
            print(f"\n💾 SUCCÈS : Données sauvegardées dans {nom_fichier}\n")
    elif profil.message_aucun_objet:
        print("Aucun objet n'a été détecté durant la session.")

    # GÉNÉRATION DES GRAPHIQUES DE VITESSE ET D'ACCÉLÉRATION
    if gestionnaire is not None:
        if archive_tableaux:
            figures_filtre(archive_tableaux, resultats["archive_filtree"], choix, gestionnaire,
                           style=profil.style_figures, fenetre_lissage=p.fenetre_lissage)
        resultats["figures"] = list(gestionnaire.fichiers)
        gestionnaire.montrer()
    elif archive_tableaux and not resultats["archive_filtree"] and profil.style_figures != STYLE_MAIN:
        seuil = "0.9" if profil.style_figures == STYLE_SAUVEGARDE else "0.5"
        print(f"Aucun objet n'a passé le filtre de ratio > {seuil}. (Graphiques propres non générés)")

    return resultats
