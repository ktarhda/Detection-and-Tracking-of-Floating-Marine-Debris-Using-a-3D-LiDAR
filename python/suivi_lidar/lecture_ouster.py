"""Lecture d'un enregistrement Ouster (fichiers .pcap et .json).

Remplace ousterFileReader, readFrame et readIMU de la Lidar Toolbox de MATLAB.
La lecture des paquets est faite par le SDK officiel d'Ouster (pip install ouster-sdk).

Conventions, choisies pour retrouver celles du code MATLAB :
- nuage organisé de taille (lignes, colonnes, 3), par exemple (128, 1024, 3), en mètres ;
- une colonne par angle d'azimut (image « destaggered ») ;
- repère « center » de ousterFileReader (valeur par défaut) : origine au centre
  optique du LiDAR, axe X vers l'angle d'encodeur 0° (côté connecteur), axe Z vers
  le haut. C'est le repère LiDAR d'Ouster. Les colonnes de la ROI (theta entre 90°
  et 270°) ont donc X négatif, comme le suppose la zone X entre -15 m et -0,4 m ;
- distance (Range) en mètres, 0 quand il n'y a pas de retour ;
- point NaN quand il n'y a pas de retour ;
- accéléromètre en m/s² et gyroscope en rad/s, lectures brutes dans les axes de l'IMU.
"""

import json

import numpy as np


def _ouvrir(constructeur, fichier, **options):
    """Ouvre une source du SDK Ouster.

    soft_id_check=True : comme ousterFileReader, on accepte les paquets même si leur
    identifiant d'initialisation diffère légèrement de celui du JSON (le SDK affiche
    alors un simple avertissement au lieu d'ignorer les paquets).
    """
    try:
        return constructeur(fichier, soft_id_check=True, **options)
    except (TypeError, ValueError):
        return constructeur(fichier, **options)


def resolution_json(fichier_json):
    """N_cols et N_rows, lus dans le JSON exactement comme dans main.m."""
    with open(fichier_json, "r", encoding="utf-8") as f:
        config = json.load(f)
    if "data_format" in config:
        format_donnees = config["data_format"]
    elif "lidar_data_format" in config:
        format_donnees = config["lidar_data_format"]
    else:
        raise ValueError("Format JSON non reconnu. Vérifiez le fichier de configuration du LiDAR.")
    return int(format_donnees["columns_per_frame"]), int(format_donnees["pixels_per_column"])


class LecteurOuster:
    """Équivalent de reader = ousterFileReader(pcapFile, jsonFile)."""

    def __init__(self, fichier_pcap, fichier_json):
        try:
            from ouster.sdk import core, pcap
        except ImportError as erreur:
            raise ImportError("Le paquet ouster-sdk est nécessaire : pip install ouster-sdk") from erreur
        self._core = core
        self._pcap = pcap
        self.fichier_pcap = str(fichier_pcap)
        self.fichier_json = str(fichier_json)

        # index=True : le fichier est indexé une fois, ce qui donne le nombre de trames
        # (reader.NumberOfFrames) et l'accès direct à une trame (readFrame(reader, i)).
        self._source = _ouvrir(pcap.PcapFrameSetSource, self.fichier_pcap, meta=[self.fichier_json], index=True)
        self.info = self._source.sensor_info[0]
        self.nombre_trames = len(self._source)
        self.lignes = self.info.format.pixels_per_column
        self.colonnes = self.info.format.columns_per_frame

        self._lut = core.XYZLut(self.info, use_extrinsics=False)       # repère capteur d'Ouster
        transformation = np.asarray(self.info.lidar_to_sensor_transform, dtype=float)
        self._R_lidar_vers_capteur = transformation[:3, :3]
        self._t_lidar_vers_capteur = transformation[:3, 3] / 1000.0    # mm -> m

    # ------------------------------------------------------------------ trames
    def _convertir(self, trame):
        """LidarFrame -> (nuage organisé dans le repère LiDAR, distances en mètres)."""
        if trame is None:
            xyz = np.full((self.lignes, self.colonnes, 3), np.nan)
            return xyz, np.zeros((self.lignes, self.colonnes))

        xyz_capteur = self._core.destagger(self.info, self._lut(trame))
        distance_mm = self._core.destagger(self.info, trame.field("RANGE"))

        # Repère capteur -> repère LiDAR : p_lidar = R' * (p_capteur - t)
        xyz = (xyz_capteur - self._t_lidar_vers_capteur) @ self._R_lidar_vers_capteur
        rho = distance_mm.astype(float) / 1000.0
        xyz[rho == 0] = np.nan
        return xyz, rho

    def lire_trame(self, i):
        """readFrame(reader, i), avec i numéroté à partir de 1 comme dans MATLAB.

        Renvoie (xyz, rho) : xyz de taille (lignes, colonnes, 3), rho de taille (lignes, colonnes).
        """
        if not 1 <= i <= self.nombre_trames:
            raise IndexError(f"Trame {i} hors du fichier (1 à {self.nombre_trames}).")
        ensemble = self._source[i - 1]
        return self._convertir(ensemble[0])

    def trames(self, premiere, derniere):
        """Parcourt les trames premiere..derniere (incluses, numérotées à partir de 1).

        Donne (i, xyz, rho) pour chaque trame, dans l'ordre. Plus rapide que des appels
        successifs à lire_trame, car le fichier est lu d'un seul trait.
        """
        premiere = max(1, int(premiere))
        derniere = min(self.nombre_trames, int(derniere))
        if derniere < premiere:
            return
        for i, ensemble in enumerate(self._source[premiere - 1:derniere], start=premiere):
            xyz, rho = self._convertir(ensemble[0])
            yield i, xyz, rho

    # ------------------------------------------------------------------ IMU
    def lire_imu(self):
        """readIMU(reader) : toutes les mesures IMU du fichier.

        Renvoie un dict :
            'accelerometre' : tableau (M, 3) en m/s² (AccelerometerReadings)
            'gyroscope'     : tableau (M, 3) en rad/s (GyroscopeReadings)
            'temps'         : tableau (M,) en secondes
        """
        source = _ouvrir(self._pcap.PcapPacketSource, self.fichier_pcap, meta=[self.fichier_json])
        accelerations, vitesses_angulaires, temps = [], [], []
        try:
            for _, paquet in source:
                if isinstance(paquet, self._core.ImuPacket):
                    accelerations.append(np.asarray(paquet.accel(), dtype=float).reshape(-1, 3))
                    vitesses_angulaires.append(np.asarray(paquet.gyro(), dtype=float).reshape(-1, 3))
                    temps.append(np.asarray(paquet.timestamp(), dtype=float).reshape(-1) * 1e-9)
        finally:
            source.close()

        if not accelerations:
            return {"accelerometre": np.empty((0, 3)), "gyroscope": np.empty((0, 3)), "temps": np.empty(0)}
        return {
            "accelerometre": np.vstack(accelerations),
            "gyroscope": np.vstack(vitesses_angulaires),
            "temps": np.concatenate(temps),
        }

    def fermer(self):
        self._source.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.fermer()
