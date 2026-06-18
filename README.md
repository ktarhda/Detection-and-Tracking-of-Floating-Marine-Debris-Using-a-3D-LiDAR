# Detection and Tracking of Floating Marine Debris Using a 3D LiDAR 

# 🌊 Experimental Validation of an Interval Particle Filter for Tracking Floating Objects

![Matlab](https://img.shields.io/badge/MATLAB-ED7D31?style=for-the-badge&logo=MathWorks&logoColor=white)
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![ROS2](https://img.shields.io/badge/ROS2-22314E?style=for-the-badge&logo=ros&logoColor=white)
![LISIC](https://img.shields.io/badge/LISIC-UR4491-blue?style=for-the-badge)

## 📌 Présentation du Projet
Ce projet est réalisé au sein du laboratoire **LISIC (UR 4491)** à l'EILCO / ULCO (Équipe EDyFI : Estimation Dynamique et Fusion d'Informations). 

L'objectif majeur est de développer, valider et implémenter une méthode robuste de **détection et de suivi d'objets flottants** (débris marins, navires, obstacles semi-immergés) en environnement aquatique. Ce système est crucial pour la sécurité maritime, la surveillance environnementale et la navigation des systèmes maritimes autonomes.

---

## 🔬 Approche Scientifique & Méthodologie

Pour faire face aux fortes incertitudes des milieux marins (vagues, reflets, conditions changeantes), le projet repose sur la fusion de capteurs et le calcul par intervalles :

1. **Module de Détection (Fusion Multi-Capteurs) :** * Fusion de données issues d'un **LiDAR 3D**, d'une **Caméra RGB** et d'une centrale **GPS-IMU**.
   * Estimation de la distance radiale et de l'orientation de l'objet sous forme d'**intervalles** afin de capturer explicitement les incertitudes de mesure.

2. **Module de Suivi (Boxed Particle Filter) :**
   * Implémentation d'un **Filtre Particulaire par Intervalles (Interval Particle Filter)**.
   * Utilisation de méthodes de Monte-Carlo basées sur des ensembles pour mettre à jour les particules sous forme de boîtes (*bounding boxes*), garantissant un suivi stable et tolérant aux ambiguïtés.

---

## 🛠️ Stack Technique & Environnements

Le projet fait le pont entre le prototypage algorithmique et la simulation robotique :

* **Langages principaux :** * **MATLAB :** Pour le prototypage mathématique, l'analyse par intervalles et la modélisation du filtre particulaire en boîte.
  * **Python :** Pour le traitement des données de vision, la manipulation des nuages de points (Point Clouds) et l'interfaçage système.
* **Architecture Mobile & Robotique :** Exploitation du framework **ROS2** (Robot Operating System) pour la gestion et la synchronisation des flux de données (*Topics*) des capteurs.
* **Environnements de Simulation & Datasets :**
  * **VRX Simulator (Virtual RobotX) :** Environnement virtuel maritime sous Gazebo pour tester les algorithmes en conditions contrôlées.
  * **Datasets Publics :** Utilisation des bases de données de référence **KITTI** et **nuScenes** pour calibrer la détection par nuages de points.
  * **Données Réelles :** Validation finale sur des scénarios réels collectés en mer/milieu aquatique par l'équipe du LISIC.

---

