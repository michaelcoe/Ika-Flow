#!/bin/bash

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
# gmshToFoam *.msh
# extrudeMesh
# createPatch -overwrite

fluentMeshToFoam -2D 1 *.msh
createPatch -dict ./system/createPatchDict_fluent -overwrite

touch mesh.foam
cd ..
