# Create the folders

rm -rf build
mkdir build
cd build

# Compile with CMake

cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8

