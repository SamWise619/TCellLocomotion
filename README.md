# 2D Bead-Spring Model of T Cell Locomotion

This repository contains a 2D bead-spring simulation code for modeling T cell migration through obstacle environments. The model represents the cell cortex and nucleus as closed bead-spring boundaries and evolves their positions under mechanical and active forces.

The code was developed for simulations of T cell locomotion, nuclear deformation, compartmentalized cytoplasmic pressure, and obstacle-mediated migration.

## Model Overview

The simulation includes the following physical mechanisms:

- cortical and nuclear stretching forces
- bending rigidity of the cortex and nucleus
- cell and nuclear area conservation
- cortex–nucleus centering force
- excluded-volume interactions between the cell membrane/cortex, nucleus, and obstacles
- rear contractile forces
- compartment boundary formation between the cortex and nucleus
- optional blebbing through oscillatory cortical spring stiffness
- optional repolarization through shifting leading-edge bead identities

## Compilation

Compile the code using `g++` with OpenMP support:

```bash
g++ -fopenmp -O3 -g -o main *.cpp
```
Running a Simulation

Run the executable with an input parameter file:

```bash
./main input.txt
```
Required Input Files

The simulation requires the following files, which are specified inside input.txt:

- positions.xyz — initial bead positions, bead types, tags, and force types
- bonds.bnd — bond connectivity and bond types
- angles.ang — angle connectivity and angle types

input.txt — simulation parameters and file paths

Example input files are included in the repository.

Key Parameters

Important parameters in input.txt include:

- dt — integration time step
- totalSimulationTime — total simulated time
- snapshotTime — time interval between saved snapshots
- K_Stretch — cortical stretching stiffness
- K_Stretch_nuc — nuclear stretching stiffness
- K_Area_CM — cortical area conservation stiffness
- K_Area_nuc — nuclear area conservation stiffness
- kCenter — cortex–nucleus centering stiffness
- kExcl_CM — excluded-volume stiffness between cell membrane/cortex and nucleus
- kExcl_Obs — excluded-volume stiffness between cell membrane/cortex and obstacles
- F_Mag_Front — magnitude of front-directed active force
- Ratio — ratio used to determine rear contractile force
- CB_On — turns compartment boundary formation on or off
- blebbing — turns oscillatory leading-edge stiffness on or off
- changeDirection — turns repolarization on or off

Hardcoded parameters not currently in input.txt:

- numNucBeads = 100beads
- numCMBeads = 200beads
- numLeadingEdgeBeads = 41beads
- d_seal = 1.3043478μm
- timeLag = 123.36953s
- updateGridTime = 1.0s
- saveInterval = 70 * snapshotStep
- calciumActivationFactor = 1.5
- repolarizationInterval = 16000.0s
- K_Stretch_min = 10000.0pN/μm
  
Output Files

The simulation writes output files named using the outputName parameter.

Typical output files include:

```bash
out-<outputName>.xyz — particle positions, bead types, and force components over time
out-<outputName>.bnd — bond information over time
```

These output files are primarily intended for visualization and analysis. They may not be directly restart-compatible with the original input file format.

Notes and Limitations

This is a research simulation code developed for a specific biophysical model. Some model choices and constants are currently hardcoded in the source code. Users should check the source carefully before modifying parameters, adding new mechanisms, or reusing the code for a different system.

Citation

If you use this code or adapt parts of it, please cite:

Alawadhi, Rutkowski, Tagay, Cartagena-Rivera, Zhovmer, Tsygankov, Vavylonis, and Tabdanov, TBD.
