# Detection and Tracking of Floating Marine Debris Using a 3D LiDAR

Master's thesis project (Master ISC, Université du Littoral Côte d'Opale), carried out at the **LISIC laboratory (UR 4491), EDyFI team**, March–September 2026.

The goal is to detect and track floating debris on the water surface with a 3D LiDAR (Ouster OS1-128) and its built-in IMU, and to compare **four state estimation filters**, probabilistic and set-membership (interval-based), on the same real data acquired on the Calais canal.

![Detection on the Calais canal](docs/detection_canal.png)

## Pipeline

```
LiDAR + IMU  ->  Tilt correction  ->  Region of interest  ->  Radial segmentation  ->  3D bounding boxes  ->  Multi-object tracking
 (Ouster)        (roll / pitch)       (water area only)       (range image, ρ)         (box centers)          (4 filters)
```

1. **Tilt correction**: roll and pitch are estimated from the IMU accelerometer, then the point cloud is rotated into the world frame at every frame.
2. **Region of interest**: only the part of the scan facing the water is kept, which reduces the number of points and removes the quay.
3. **Radial segmentation**: the range matrix ρ (128 × 512) is processed as a depth image. A jump between neighbouring cells larger than an adaptive threshold τ = max(ε, 1.5 μ) marks the edge of an object. Regions are then closed (3×3) and labelled (8-connectivity).
4. **Bounding boxes**: one axis-aligned 3D box per object; its center is the measurement sent to the tracker.
5. **Multi-object tracking**: constant-acceleration model (9 states: position, velocity, acceleration), Mahalanobis gating (d_M < 3.5) and nearest-neighbour association.

## The four filters

| Filter | Family | Uncertainty representation | File |
| --- | --- | --- | --- |
| Kalman filter | Probabilistic (Gaussian) | Mean + covariance | `kalman_classique.m` |
| Particle filter | Probabilistic (Monte-Carlo) | Weighted particles + resampling | `filtre_particule.m` |
| UBIKF (Upper Bound Interval Kalman Filter) | Set-membership | State box + upper bound of the covariance | `kalman_ensembliste.m` |
| Box Particle Filter | Set-membership (Monte-Carlo) | Weighted boxes + subdivision resampling | `filtre_particule_boite.m` |

## Results (Calais canal, 3 floating objects)

<table>
  <tr>
    <td align="center"><img src="docs/canal_calais.png" height="320"><br><em>Three floating objects on the Calais canal</em></td>
    <td align="center"><img src="docs/suivi_objets.png" height="320"><br><em>Detection and tracking: each object keeps its box and its ID</em></td>
  </tr>
</table>

| Filter | Mean number of tracks | Mean cardinality error | Frames with exact count |
| --- | --- | --- | --- |
| Kalman | 3.04 | 0.10 | 92.1 % |
| Particle filter | 2.94 | 0.06 | 98.1 % |
| UBIKF | 2.94 | 0.06 | 98.1 % |
| Box Particle Filter | 2.96 | 0.07 | 96.8 % |

All four filters estimate drift speeds between 0.12 and 0.25 m/s, consistent with the current measured on site (about 0.15 m/s).

![3D trajectories of object ID 1 for the four filters](docs/trajectoires.png)

## Repository structure

```
├── main.m                      # Kalman filter or UBIKF (choose at the top of the file)
├── Detect_particule.m          # Particle filter or Box Particle Filter (choose at the top of the file)
├── correct_inclinaison.m       # Tilt correction with the IMU
├── segmentation_rho_marin.m    # Radial segmentation and bounding boxes
├── kalman_classique.m          # Kalman filter
├── kalman_ensembliste.m        # UBIKF
├── filtre_particule.m          # Particle filter
├── filtre_particule_boite.m    # Box Particle Filter
├── sauvegarde.m                # Saves the results of a run (.mat)
├── plot_figure.m               # Compares the four filters from the saved .mat files
└── docs/                       # Figures, report and slides
```

## Requirements

- MATLAB (developed with R2024b)
- Lidar Toolbox (reading Ouster `.pcap` files with `ousterFileReader`)
- Image Processing Toolbox (morphological closing, connected-component labelling)
- [INTLAB](https://www.tuhh.de/ti3/rump/intlab/) for interval computations (UBIKF and Box Particle Filter), not included in this repository

## How to run

1. Install INTLAB and add it to the MATLAB path.
2. Put the acquisition files (`3objets.pcap` and `3objets.json`) in the project folder. The raw LiDAR recordings are not included in this repository because of their size.
3. Run `main.m` for the Kalman filter or the UBIKF, or `Detect_particule.m` for the particle filter or the Box Particle Filter. The filter is selected by the `CHOIX_FILTRE` variable at the top of each script.
4. Run `sauvegarde.m` to save the results, then `plot_figure.m` to compare the four filters.

## References

- T. A. Tran, C. Jauberthie, L. Travé-Massuyès, Q. H. Lu, *An Interval Kalman Filter enhanced by lowering the covariance matrix upper bound*, International Journal of Applied Mathematics and Computer Science, 31(2), 2021.
- J. Xiong, C. Jauberthie, L. Travé-Massuyès, F. Le Gall, *Fault Detection using Interval Kalman Filtering enhanced by Constraint Propagation*, IEEE CDC, 2013.
- F. Abdallah, A. Gning, P. Bonnifait, *Box particle filtering for nonlinear state estimation using interval analysis*, Automatica, 44(3), 2008.
- S. M. Rump, *INTLAB – INTerval LABoratory*, Developments in Reliable Computing, 1999.

## Author

**Khalil Tarhda**, supervised by Régis Lherbier and Mohamed Fnadi (LISIC, Université du Littoral Côte d'Opale).
