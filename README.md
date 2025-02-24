# Ika Flow

**Ika** is a Māori word meaning any creature that swims in fresh or salt water including marine mammals such as whales.

This is a repository for code to estimate power and thrust of an undulating fish and related simulations.

## Compiling the solver

Copy the folder OpenFOAM -> dynamicMesh to your working directory. Run wclean and wmake to make compile the solver. Note that the solver is working with OpenFoam v2206, currently.

## TestCase

A test case is supplied in the folder OpenFOAM -> testCase. There are two folders here: mesh, and st0_40. The mesh folder holds the mesh and background mesh. The st0_40 has carangiform motion at a Strouhal number of 40.

First, you need to make sure all the files are executable. You can navigate to the top level folder **testCase** and run the following command:

```console
find . -wholename "**/*.sh" -exec chmod +x {} \;
```

Then you can create the mesh by navigating to the OpenFOAM -> testCase -> mesh folder and running

```console
./run_mesh.sh
```

The important part here is that the **topSetDict_movingZone** is executed. This is what sets the overset mesh as the moving zone and lets the solver know you want to move this part of the mesh.

To run the case, you navigate to OpenFOAM -> testCase -> st0_40 and use the command

```console
./run_all.sh
```

Important here is that the ***controlDict** has libs entries for both "liboverset.so" and the solver entry (which I have named "libfishBodyMotion.so"). Another important point here is that the dynamicMeshDict is set up with the moving zone and the proper coefficients to define the movement you want.

## File Structure

The folders are structured as follows:

### GMSH

The GMSH folder has scrips and **.geo** files for making an airfoil in gmsh.

### Lighthill

The Lighthill folder contains ipython notebooks to model forces and energy via Lighthill's theorem. This was just a play around and shouldn't be taken too seriously.

### Motion Modeling

This folder was used to model motions, forces, and energy from images of fish.

### OpenFOAM

This folder holds all the files for the solver and a test case. Furthermore, all the post-processing scripts for after the cases are run.

### Reference_data_modeling

This folder contains the script to compare the results to that of Yu et al.