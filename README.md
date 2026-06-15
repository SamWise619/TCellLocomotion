# TCellLocomotion

2D bead-spring code that simulates cell migration in obstacle environments by evolving the cell and nucleus boundaries according to forces of stretching, bending, centering, contraction, and area conservation under the effect of a nuclear piston. This code was developed for the simulations in Alawadhi, Rutkowski, Tagay, Cartagena-Rivera, Zhovmer, Tysgankov, Vavylonis, and Tabdanov, TBD.

Compilation using g++: g++ -fopenmp -O3 -g -o main *.cpp

Input parameter file (input.txt): This file lists various parameters used by the simulation. Command to run code: ./main input.txt

Must supply positions.xyz, bonds.bnd, and angles.ang files for the input position, bond, and angle files, respectively, into input.txt (see example files).

# Revision: 2D Bead-Spring Model of T Cell Locomotion

This repository contains a 2D bead-spring simulation code for modeling T cell migration through obstacle environments. The model represents the cell cortex and nucleus as closed bead-spring boundaries and evolves their positions under mechanical and active forces.

The code was developed for simulations of T cell locomotion, nuclear deformation, compartmentalized cytoplasmic pressure, and obstacle-mediated migration.

## Model Overview

The simulation includes the following physical mechanisms:

- cortical and nuclear stretching forces
- bending rigidity of the cortex and nucleus
- cell and nuclear area conservation
- cortex–nucleus centering force
- excluded-volume interactions between the cell, nucleus, and obstacles
- rear contractile forces
- leading-edge force generation
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

positions.xyz — initial bead positions, bead types, tags, and force types
bonds.bnd — bond connectivity and bond types
angles.ang — angle connectivity and angle types

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

Output Files

The simulation writes output files named using the outputName parameter.

Typical output files include:

- out-<outputName>.xyz — particle positions, bead types, and force components over time
- out-<outputName>.bnd — bond information over time

These output files are primarily intended for visualization and analysis. They may not be directly restart-compatible with the original input file format.

Notes and Limitations

This is a research simulation code developed for a specific biophysical model. Some model choices and constants are currently hardcoded in the source code. Users should check the source carefully before modifying parameters, adding new mechanisms, or reusing the code for a different system.

Citation

If you use this code or adapt parts of it, please cite:

Alawadhi, Rutkowski, Tagay, Cartagena-Rivera, Zhovmer, Tsygankov, Vavylonis, and Tabdanov, TBD.
