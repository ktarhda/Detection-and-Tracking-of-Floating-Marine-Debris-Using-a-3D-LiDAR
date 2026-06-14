# Detection and Tracking of Floating Marine Debris Using a 3D LiDAR 

# 🌊 Experimental Validation of an Interval Particle Filter for Tracking Floating Objects

![Matlab](https://img.shields.io/badge/MATLAB-ED7D31?style=for-the-badge&logo=MathWorks&logoColor=white)
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![ROS2](https://img.shields.io/badge/ROS2-22314E?style=for-the-badge&logo=ros&logoColor=white)
![LISIC](https://img.shields.io/badge/LISIC-UR4491-blue?style=for-the-badge)

## 📌 Présentation du Projet
[cite_start]Ce projet est réalisé au sein du laboratoire **LISIC (UR 4491)** à l'EILCO / ULCO (Équipe EDyFI : Estimation Dynamique et Fusion d'Informations)[cite: 3, 655]. 

[cite_start]L'objectif majeur est de développer, valider et implémenter une méthode robuste de **détection et de suivi d'objets flottants** (débris marins, navires, obstacles semi-immergés) en environnement aquatique[cite: 630]. [cite_start]Ce système est crucial pour la sécurité maritime, la surveillance environnementale et la navigation des systèmes maritimes autonomes[cite: 630].

---

## 🔬 Approche Scientifique & Méthodologie

Pour faire face aux fortes incertitudes des milieux marins (vagues, reflets, conditions changeantes), le projet repose sur la fusion de capteurs et le calcul par intervalles :

1. [cite_start]**Module de Détection (Fusion Multi-Capteurs) :** * Fusion de données issues d'un **LiDAR 3D**, d'une **Caméra RGB** et d'une centrale **GPS-IMU**[cite: 631, 636].
   * [cite_start]Estimation de la distance radiale et de l'orientation de l'objet sous forme d'**intervalles** afin de capturer explicitement les incertitudes de mesure[cite: 636].

2. **Module de Suivi (Boxed Particle Filter) :**
   * [cite_start]Implémentation d'un **Filtre Particulaire par Intervalles (Interval Particle Filter)**[cite: 632].
   * [cite_start]Utilisation de méthodes de Monte-Carlo basées sur des ensembles pour mettre à jour les particules sous forme de boîtes (*bounding boxes*), garantissant un suivi stable et tolérant aux ambiguïtés[cite: 634, 638].

---

## 🛠️ Stack Technique & Environnements

Le projet fait le pont entre le prototypage algorithmique et la simulation robotique :

* [cite_start]**Langages principaux :** * **MATLAB :** Pour le prototypage mathématique, l'analyse par intervalles et la modélisation du filtre particulaire en boîte.
  * [cite_start]**Python :** Pour le traitement des données de vision, la manipulation des nuages de points (Point Clouds) et l'interfaçage système.
* [cite_start]**Architecture Mobile & Robotique :** Exploitation du framework **ROS2** (Robot Operating System) pour la gestion et la synchronisation des flux de données (*Topics*) des capteurs[cite: 657].
* **Environnements de Simulation & Datasets :**
  * [cite_start]**VRX Simulator (Virtual RobotX) :** Environnement virtuel maritime sous Gazebo pour tester les algorithmes en conditions contrôlées[cite: 645, 648].
  * [cite_start]**Datasets Publics :** Utilisation des bases de données de référence **KITTI** et **nuScenes** pour calibrer la détection par nuages de points[cite: 650].
  * [cite_start]**Données Réelles :** Validation finale sur des scénarios réels collectés en mer/milieu aquatique par l'équipe du LISIC[cite: 651].

---

