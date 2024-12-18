#!/bin/bash -e

foamCleanTutorials

cd airfoil
sh run_mesh.sh
cd ..

cd backGround
blockMesh
#snappyHexMesh -overwrite

topoSet -dict system/topoSetDictR1

refineMesh -dict system/refineMeshDict1 -overwrite

topoSet -dict system/topoSetDictR2

refineMesh -dict system/refineMeshDict2 -overwrite

mergeMeshes . ../airfoil -overwrite

topoSet
topoSet -dict system/topoSetDict_movingZone

checkMesh -allTopology -allGeometry -writeSets vtk |  tee log.checkMesh

touch mesh.foam

cd ..
