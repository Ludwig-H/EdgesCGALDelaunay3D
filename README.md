# EdgesCGALDelaunay3D

Extraction rapide du **1‑squelette** de la triangulation de **Delaunay 3D** avec **CGAL**.

- **Entrée**: `.xyz` (une ligne: `x y z`)
- **Sortie**: fichier texte, une ligne par arête finie: `i j` (indices 0‑based)

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

```bash
./build/EdgesCGALDelaunay3D data/example.xyz out.edges
# Parallélisme (si TBB + tbbmalloc détectés) :
CGAL_NTHREADS=$(nproc) ./build/EdgesCGALDelaunay3D data/example.xyz out.edges
```

## Notes
- Lecture/écriture bufferisées, insertion **par lot** via `(Point, index)` pour vitesse et faible overhead.
- Mode parallèle activé uniquement si `tbb` et `tbbmalloc` sont trouvés pour éviter les erreurs de linkage.
