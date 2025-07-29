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

Build version
```
cmake -DCGAL_DIR=/home/peter/Projects/libraries/cgal/build -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On ..
make
```

Release version
```
cmake -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/cgal/install;/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..
make
```


# Test

make && ./fv99 -f ../data/three-sheet-toy.txt
make && ./fv99 -f ~/Projects/data/reeb-space-test-data/data.vtu
make && ./fv99 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-50-tets-320.vtu -e 1e-3
make && ./fv99 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-20-300.vtu
