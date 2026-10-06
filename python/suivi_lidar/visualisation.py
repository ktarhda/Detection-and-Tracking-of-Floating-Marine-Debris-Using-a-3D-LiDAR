"""Affichage 3D en direct, en remplacement de pcplayer (blocs COLORATION et AFFICHAGE de sauvegarde.m).

Le nuage est coloré par la distance (palette jet, 40 m maximum), les points situés
dans une boîte de détection sont en bleu, les boîtes de détection en vert, et chaque
piste confirmée est dessinée avec sa trajectoire estimée, sa trajectoire observée
(jaune), sa position actuelle et son texte « ID | vitesse ». La boîte UBIKF est en
blanc, les boîtes du BPF en orange et les particules en gris clair.

L'affichage utilise Open3D (pip install open3d). Fermer la fenêtre arrête la boucle
principale, comme la fermeture du pcplayer dans MATLAB.
"""

import numpy as np

from .parametres import BOX_PARTICULAIRE, CLASSIQUE, ENSEMBLISTE, PARTICULAIRE


def jet_matlab(m=256):
    """jet(m) de MATLAB (même définition que la fonction jet.m)."""
    n = int(np.ceil(m / 4))
    u = np.concatenate([np.arange(1, n + 1) / n, np.ones(n - 1), np.arange(n, 0, -1) / n])
    g = int(np.ceil(n / 2)) - (1 if m % 4 == 1 else 0) + np.arange(1, len(u) + 1)
    r = g + n
    b = g - n
    g, r, b = g[g <= m], r[r <= m], b[b >= 1]
    J = np.zeros((m, 3))
    J[r - 1, 0] = u[:len(r)]
    J[g - 1, 1] = u[:len(g)]
    J[b - 1, 2] = u[len(u) - len(b):]
    return J


_PALETTE = jet_matlab(256)


def couleurs_nuage(rho_brut, xyz_aligne, boites, distance_max_couleur=40.0):
    """Bloc COLORATION : couleurs uint8 (lignes, colonnes, 3).

    Couleur jet selon la distance, puis bleu pour les points dans une boîte de détection.
    """
    rho_norm = np.minimum(rho_brut / distance_max_couleur, 1.0)
    indices = np.floor(rho_norm * 255 + 0.5).astype(int)          # round(rho_norm*255) + 1 en MATLAB
    couleurs = np.floor(_PALETTE[indices] * 255 + 0.5).astype(np.uint8)

    X, Y, Z = xyz_aligne[..., 0], xyz_aligne[..., 1], xyz_aligne[..., 2]
    with np.errstate(invalid="ignore"):
        for box in np.asarray(boites).reshape(-1, 9):
            masque = ((X >= box[0] - box[3] / 2) & (X <= box[0] + box[3] / 2)
                      & (Y >= box[1] - box[4] / 2) & (Y <= box[1] + box[4] / 2)
                      & (Z >= box[2] - box[5] / 2) & (Z <= box[2] + box[5] / 2))
            couleurs[masque] = (0, 0, 255)
    return couleurs


def points_affichables(xyz_aligne, couleurs):
    """select(ptCloud_visu, idx_valides) : points dans [-50, 50] x [-50, 50] x [-5, 5]."""
    xyz = xyz_aligne.reshape(-1, 3)
    c = couleurs.reshape(-1, 3)
    with np.errstate(invalid="ignore"):
        valides = ((xyz[:, 0] >= -50) & (xyz[:, 0] <= 50) & (xyz[:, 1] >= -50) & (xyz[:, 1] <= 50)
                   & (xyz[:, 2] >= -5) & (xyz[:, 2] <= 5))
    return xyz[valides], c[valides]


# Couleurs des trajectoires et de la position actuelle (sauvegarde.m)
COULEUR_FILTRE = {
    CLASSIQUE: (1.0, 0.0, 0.0),          # 'r'
    ENSEMBLISTE: (0.0, 1.0, 0.0),        # 'g'
    PARTICULAIRE: (0.0, 1.0, 1.0),       # 'c'
    BOX_PARTICULAIRE: (1.0, 0.5, 0.0),   # orange
}
JAUNE = (1.0, 1.0, 0.0)
VERT = (0.0, 1.0, 0.0)
BLANC = (1.0, 1.0, 1.0)


class Visualiseur:
    """Fenêtre 3D mise à jour à chaque trame (équivalent de viewer = pcplayer(...))."""

    def __init__(self, choix_filtre, trace_observe=True, largeur=1280, hauteur=800, taille_texte=0.012):
        try:
            import open3d as o3d
        except ImportError as erreur:
            raise ImportError("L'affichage 3D demande Open3D : pip install open3d") from erreur
        self.o3d = o3d
        self.choix_filtre = choix_filtre
        self.trace_observe = trace_observe      # trajectoire observée en jaune (pas dans Detect_particule.m)
        self.taille_texte = taille_texte
        self.vis = o3d.visualization.Visualizer()
        if not self.vis.create_window(window_name=f"Détection Objets Aquatiques : {choix_filtre}",
                                      width=largeur, height=hauteur):
            raise RuntimeError("Impossible d'ouvrir la fenêtre 3D (pas d'écran disponible ?).")
        options = self.vis.get_render_option()
        options.background_color = np.array([0.0, 0.0, 0.0])
        options.point_size = 1.5
        options.line_width = 2.0
        options.mesh_show_back_face = True          # le texte reste visible des deux côtés
        self.nuage = o3d.geometry.PointCloud()
        self.geometries = []
        self.premiere_trame = True
        self.ouvert = True

    # ------------------------------------------------------------------ outils de dessin
    def _ajouter(self, geometrie):
        self.vis.add_geometry(geometrie, reset_bounding_box=False)
        self.geometries.append(geometrie)

    def _boite(self, centre, dimensions, couleur):
        o3d = self.o3d
        centre = np.asarray(centre, dtype=float)
        demi = np.maximum(np.asarray(dimensions, dtype=float), 1e-3) / 2
        if not np.all(np.isfinite(centre)) or not np.all(np.isfinite(demi)):
            return
        aabb = o3d.geometry.AxisAlignedBoundingBox(centre - demi, centre + demi)
        lignes = o3d.geometry.LineSet.create_from_axis_aligned_bounding_box(aabb)
        lignes.paint_uniform_color(couleur)
        self._ajouter(lignes)

    def _ligne(self, points, couleur):
        """Trace une polyligne ; les NaN coupent la ligne, comme plot3 dans MATLAB."""
        o3d = self.o3d
        points = np.asarray(points, dtype=float)
        if len(points) < 2:
            return
        valides = np.all(np.isfinite(points), axis=1)
        segments = np.flatnonzero(valides[:-1] & valides[1:])
        if segments.size == 0:
            return
        lignes = o3d.geometry.LineSet(o3d.utility.Vector3dVector(points),
                                      o3d.utility.Vector2iVector(np.column_stack([segments, segments + 1])))
        lignes.paint_uniform_color(couleur)
        self._ajouter(lignes)

    def _points(self, points, couleur):
        o3d = self.o3d
        points = np.asarray(points, dtype=float)
        points = points[np.all(np.isfinite(points), axis=1)]
        if len(points) == 0:
            return
        nuage = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(points))
        nuage.paint_uniform_color(couleur)
        self._ajouter(nuage)

    def _sphere(self, centre, couleur, rayon=0.08):
        if not np.all(np.isfinite(centre)):
            return
        sphere = self.o3d.geometry.TriangleMesh.create_sphere(radius=rayon, resolution=8)
        sphere.translate(np.asarray(centre, dtype=float))
        sphere.paint_uniform_color(couleur)
        self._ajouter(sphere)

    def _texte(self, texte, position):
        """Texte blanc 0,8 m au-dessus de l'objet, lisible depuis le LiDAR (regard vers -X)."""
        if not np.all(np.isfinite(position)):
            return
        try:
            maillage = self.o3d.t.geometry.TriangleMesh.create_text(texte, depth=0.0).to_legacy()
        except Exception:
            return
        maillage.scale(self.taille_texte, center=(0.0, 0.0, 0.0))
        # Axes du texte : x -> +Y du monde, y -> +Z du monde
        maillage.rotate(np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), center=(0.0, 0.0, 0.0))
        maillage.translate(np.asarray(position, dtype=float) + np.array([0.0, 0.0, 0.8]))
        maillage.paint_uniform_color(BLANC)
        self._ajouter(maillage)

    # ------------------------------------------------------------------ mise à jour
    def est_ouvert(self):
        """isOpen(viewer)."""
        return self.ouvert

    def afficher(self, xyz, couleurs_uint8, boites, pistes_affichees):
        """view(viewer, ptCloud_visu) puis dessin des boîtes et des pistes.

        xyz, couleurs_uint8 : points déjà sélectionnés (n, 3)
        boites              : boîtes de détection (n, 9)
        pistes_affichees    : liste de (piste, etat) pour les pistes confirmées
        """
        o3d = self.o3d
        for geometrie in self.geometries:
            self.vis.remove_geometry(geometrie, reset_bounding_box=False)
        self.geometries = []

        self.nuage.points = o3d.utility.Vector3dVector(np.asarray(xyz, dtype=float))
        self.nuage.colors = o3d.utility.Vector3dVector(np.asarray(couleurs_uint8, dtype=float) / 255.0)
        if self.premiere_trame:
            self.vis.add_geometry(self.nuage)
            vue = self.vis.get_view_control()
            vue.set_lookat([-7.0, 0.0, -1.0])         # centre de la zone d'eau
            vue.set_front([0.8, 0.0, 0.6])            # caméra derrière et au-dessus du LiDAR
            vue.set_up([0.0, 0.0, 1.0])
            # Le zoom d'Open3D dépend de la taille du nuage : on vise une vue d'environ 15 m de large.
            diagonale = float(np.linalg.norm(self.nuage.get_max_bound() - self.nuage.get_min_bound()))
            vue.set_zoom(float(np.clip(6.0 / max(diagonale, 1.0), 0.02, 1.0)))
            self.premiere_trame = False
        else:
            self.vis.update_geometry(self.nuage)

        # Boîtes de détection (showShape('cuboid', boites, 'Color', [0 1 0]))
        for box in np.asarray(boites).reshape(-1, 9):
            self._boite(box[0:3], box[3:6], VERT)

        for piste, etat in pistes_affichees:
            couleur = COULEUR_FILTRE[self.choix_filtre]
            if self.choix_filtre == ENSEMBLISTE:
                boite = etat["boite_position"]
                self._boite(boite.inf + (boite.sup - boite.inf) / 2, boite.sup - boite.inf, BLANC)
            elif self.choix_filtre == PARTICULAIRE:
                self._points(piste["particules"][0:3].T, (0.5, 0.5, 0.5))
            elif self.choix_filtre == BOX_PARTICULAIRE:
                boites_piste = etat["boites"]
                for b in range(boites_piste.shape[1]):
                    bas, haut = boites_piste.inf[0:3, b], boites_piste.sup[0:3, b]
                    dimensions = np.maximum(0.01, haut - bas)
                    opacite = min(1.0, max(0.1, float(etat["poids"][b]) * 2))
                    self._boite(bas + dimensions / 2, dimensions, tuple(opacite * np.array(couleur)))

            historique = piste["historique"]
            if len(historique) > 1:
                if self.choix_filtre == CLASSIQUE:
                    self._points(historique, couleur)       # 'r*'
                else:
                    self._ligne(historique, couleur)
            if self.trace_observe and len(piste["historique_observe"]) > 1:
                self._ligne(piste["historique_observe"], JAUNE)
            self._sphere(etat["pos_actuelle"], couleur)
            self._texte(etat["texte_affichage"], etat["pos_actuelle"])

        self.ouvert = bool(self.vis.poll_events())
        self.vis.update_renderer()

    def capture(self, chemin):
        """Enregistre une image de la fenêtre (utile pour un rapport)."""
        self.vis.capture_screen_image(str(chemin), do_render=True)

    def fermer(self):
        try:
            self.vis.destroy_window()
        except Exception:
            pass
        self.ouvert = False
