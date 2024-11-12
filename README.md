# Description

This software is made of MATLAB code accompanying the manuscript

> *Robust undulatory locomotion via neuromechanical adjustments\\ in a dissipative medium* by K. Ishimoto, C. Moreau, and J. Herault
>
accepted in Journal of Royal Society Interface as of Nov. 12, 2024.

# Use

Minimal code for exploring the dynamics is contained in the "main_CPG.m" file.
Dependent functions are in the "helpers" folder.

## "helpers" folder

> "coordinates.filament" : computes the coordinates of every hinge from the position and angle data.
> "cpg_sym": computes the right-hand side of the dynamical system to be solved.
>  "matrix3Nparameters","matrixTorque" : helpers for computing the system that needs to be inverted (the hydrodynamic matrix A)

