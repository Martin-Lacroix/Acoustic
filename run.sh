# Path to the input files

export MESH=${PWD}/examples/2DSound.msh
export CONFIG=${PWD}/examples/config.conf

# Clean output folder

rm -rf workspace
mkdir workspace
cd workspace

# Run the code

export OMP_NUM_THREADS=8
../build/Acoustic ${MESH} ${CONFIG} 2>&1 | tee workspace.txt
