

import numpy as np

# Translations du fichier JSON, écrites exactement comme dans correct_inclinaison.m.
# Remarque : les valeurs de t_imu_sensor sont déjà en mètres, la division par 1000
# les rend donc presque nulles. L'effet sur le nuage est inférieur au millimètre ;
# on garde le calcul MATLAB pour obtenir exactement les mêmes valeurs.
T_IMU_SENSOR = np.array([0.006253, -0.011775, 0.007645]) / 1000
T_LIDAR_SENSOR = np.array([0.0, 0.0, 36.18]) / 1000


def Rx(roll):
    c, s = np.cos(roll), np.sin(roll)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, c, -s],
                     [0.0, s, c]])


def Ry(pitch):
    c, s = np.cos(pitch), np.sin(pitch)
    return np.array([[c, 0.0, s],
                     [0.0, 1.0, 0.0],
                     [-s, 0.0, c]])


def moyenne_accelerometre(imu):
    """Renvoie [ax, ay, az].

    imu : dict renvoyé par LecteurOuster.lire_imu() (clé 'accelerometre', tableau N x 3),
          tableau N x 3, ou vecteur [ax, ay, az] déjà moyenné (cas « temps réel » de MATLAB).
    """
    if isinstance(imu, dict):
        imu = imu["accelerometre"]
    acc = np.asarray(imu, dtype=float)
    if acc.ndim == 2:
        return acc.mean(axis=0)          # mean(accel_data(:,1)), ...
    return acc[:3]


def roll_pitch(imu):
    """Roll et pitch (radians), comme dans correct_inclinaison.m."""
    ax, ay, az = moyenne_accelerometre(imu)
    roll = np.arctan2(ay, az)
    pitch = np.arctan2(-ax, np.sqrt(ay ** 2 + az ** 2))
    return roll, pitch


def matrice_correction(imu):
    """R_correction = (Ry * Rx)'."""
    roll, pitch = roll_pitch(imu)
    return (Ry(pitch) @ Rx(roll)).T


def correct_inclinaison(points, imu):
    """Remet le nuage à l'horizontale.

    points : tableau (..., 3) en mètres (lignes x colonnes x 3 pour un nuage organisé).
             Les NaN (pas de retour) restent NaN.
    imu    : voir moyenne_accelerometre.
    Renvoie un tableau de même forme (ptCloud_corrige).
    """
    R_correction = matrice_correction(imu)

    # Décalage entre le centre LiDAR et l'IMU
    decalage_lidar_imu = T_IMU_SENSOR - T_LIDAR_SENSOR

    p = np.asarray(points, dtype=float)
    forme = p.shape
    p = p.reshape(-1, 3)
    p = p - decalage_lidar_imu           # T1 : déplacer le nuage au centre de l'IMU
    p = p @ R_correction.T               # T2 : rotation de correction (p' = R * p)
    p = p + decalage_lidar_imu           # T3 : ramener au centre du LiDAR
    return p.reshape(forme)
