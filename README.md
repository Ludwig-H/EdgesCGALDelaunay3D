# EdgesCGALDelaunay3D

Extraction rapide du **1‑squelette** de la triangulation de **Delaunay 3D** avec **CGAL**.

- **Entrée**: `.npy` binaire `(N,3)` (`float32` ou `float64`)
- **Sortie**: `.npy` binaire `(M,2)` (`uint32` indices 0‑based)

## Dépendances (Ubuntu 22.04/24.04)

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libcgal-dev libtbb-dev libtbbmalloc2 libgmp-dev libmpfr-dev
```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Utilisation
Fournissez un fichier `.npy` contenant un tableau `(N,3)` de points.

```bash
./build/EdgesCGALDelaunay3D points.npy edges.npy
# Parallélisme (si TBB + tbbmalloc détectés) :
CGAL_NTHREADS=$(nproc) ./build/EdgesCGALDelaunay3D points.npy edges.npy
```

Pour générer un nuage de points d'exemple :

```bash
python - <<'PY'
import numpy as np
np.save('points.npy', np.random.rand(100,3).astype(np.float32))
PY
```

## Notes
- Lecture/écriture bufferisées, insertion **par lot** via `(Point, index)` pour vitesse et faible overhead.
- Mode parallèle activé uniquement si `tbb` et `tbbmalloc` sont trouvés pour éviter les erreurs de linkage.
