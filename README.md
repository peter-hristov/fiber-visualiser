# Compile

Here's how to compile the code with cmake, using CGAL version 6.0.1 and VTK 9.4.1

```
cmake -DCGAL_DIR=/home/peter/Projects/libraries/cgal/build -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On ..

make


```



# Test

make && ./fv99 ../data/vertex-neighbourhood-7-vertices.vtu
make && ./fv99 ../data/three-sheet-toy.txt
