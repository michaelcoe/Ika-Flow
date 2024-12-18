#!/bin/bash -e
#SBATCH --ntasks 8
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 150:00:00
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 94:00:00
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 94:00:00
#SBATCH --export NONE

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 94:00:00
#SBATCH --export NONE

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 94:00:00
#SBATCH --export NONE

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --mem-per-cpu 1500
#SBATCH --time 94:00:00
#SBATCH --export NONE

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 85:00:00
#SBATCH --export NONE

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 85:00:00
#SBATCH --export NONE

module load OpenFOAM/v2106-gimkl-2020a && source $FOAM_BASH

processors=32 
 
cd backGround

rm -r 0
cp -r 0_org 0

setFields | tee log.setFields

setFields -dict system/setFieldsDict_init 

decomposePar -force | tee log.decomposePar

mpirun -np $processors renumberMesh -overwrite -parallel | tee log.renumberMesh
mpirun -np $processors overPimpleDyMFoam -parallel | tee log.solver

rm -r processor*

#renumberMesh -overwrite | tee log.renumberMesh
#overPimpleDyMFoam | tee log.solver
echo "$(pwd)" >> ~/OpenFOAM/rccuser-v2106/run/carangiform_amp_2_00_k_1_00/completed.txt
