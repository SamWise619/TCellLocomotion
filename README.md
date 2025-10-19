# TCellLocomotion

2D bead-spring code that simulates cell migration in obstacle environments by evolving the cell and nucleus boundaries according to forces of stretching, bending, centering, contraction, and area conservation under the effect of a nuclear piston. This code was developed for the simulations in Alawadhi, Rutkowski, Tagay, Cartagena-Rivera, Zhovmer, Tysgankov, Vavylonis, and Tabdanov, TBD.

Compilation using g++: g++ -fopenmp -O3 -g -o main *.cpp

Input parameter file (input.txt): This file lists various parameters used by the simulation. Command to run code: ./main input.txt

Must supply positions.xyz, bonds.bnd, and angles.ang files for the input position, bond, and angle files, respectively, into input.txt (see example files).
