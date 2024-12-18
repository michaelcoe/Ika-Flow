#!/bin/bash -e

decomposePar -force
mpirun -np 4 renumberMesh -overwrite -parallel
mpirun -np 4 overPimpleDyMFoam -parallel

rm -r processor*