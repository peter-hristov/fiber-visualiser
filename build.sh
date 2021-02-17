# Store the project directory
DIRECTORY=$(cd `dirname $0` && pwd)

# Add submodules
cd $DIRECTORY
git submodule update --init --recursive

# Make netcdf build and install directories
cd $DIRECTORY/external/build
rm -rf netcdf-c
mkdir -p netcdf-c/build netcdf/install

# Build and install netcdf
cd $DIRECTORY/external/build/netcdf-c/build
cmake -D CMAKE_INSTALL_PREFIX="$DIRECTORY/external/build/netcdf-c/install" "$DIRECTORY/external/netcdf-c"
make
make install

## Build tv9k
cd $DIRECTORY/build
rm -rf *
cmake -D CMAKE_PREFIX_PATH="$DIRECTORY/external/build/netcdf-c/install" ..
make
