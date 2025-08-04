# Arrange and Traverse Algorithm for Exact Reeb Space Computation
An open-source library and visualisation tool for exact Reeb space computation. This is based on on the paper Arrange and Traverse Algorithm for Computation of Reeb Spaces of Piecewise Linear Maps, presented at SGP 2025. The manuscript is published by CGF - https://diglib.eg.org/items/5c5acfeb-1389-4ea5-a036-ce9282ab4bab. Please cite this paper as.

```
@article{hristov2025,
    journal = {Computer Graphics Forum},
    title = {{Arrange and Traverse Algorithm for Computation of Reeb Spaces of Piecewise Linear Maps}},
    author = {Hristov, Petar and Sakurai, Daisuke and Carr, Hamish and Hotz, Ingrid and Masood, Talha Bin},
    year = {2025},
    publisher = {The Eurographics Association and John Wiley & Sons Ltd.},
    ISSN = {1467-8659},
    DOI = {10.1111/cgf.70206}
}
        
```

This code is distributed under the terms of the GNU General Public License version 3 (GPLv3). You are free to use, modify, and redistribute it, provided that any derivative works also remain under the GPLv3 license. For more details see ./LICENSE.


# Installing Dependencies

## CGAL

```
git clone --recursive https://github.com/CGAL/cgal cgal
cd ./cgal
git checkout v6.0.1 

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/librarires/cgal/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install
```


## VTK

```
git clone --recursive https://gitlab.kitware.com/vtk/vtk.git VTK-9.4.1
cd ./VTK-9.4.1
git checkout v9.4.1 

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/librarires/VTK-9.4.1/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install
```


# Compiling

## Build version
```
cmake -DCGAL_DIR=/home/peter/Projects/libraries/cgal/build -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On ..
make
```

## Release version
```
cmake -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/cgal/install;/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..
make
```


# Testing

```
cd ./build
make && ./fv99 -f ../data/three-sheet-toy.txt
make && ./fv99 -f ../data/vertex-neighbourhood-7-vertices.vtu
make && ./fv99 -f ../data/torus-factor-50-tets-320.vtu -e 1e-3
make && ./fv99 -f ../data//downsample-20-300.vtu
```
