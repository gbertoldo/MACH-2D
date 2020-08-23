# MACH-2D 5.10.2
2D flow field simulation.

Solves external and internal compressible flow.

Bugs fixed relatively to version 5.10.1:
- Temperature and pressure were not initialized in fictitious corners
neither on boundary faces.
- Removed the conditional that forbade application of total enthalpy
conservation to calculate the temperature in the internal flow.
- Fixed the application of boundary conditions for Navier-Stokes equation in the internal flow.

Improvements relatively to version 5.10.1:
- For Euler model with constant thermophysical properties, T is
calculated from conservation of total enthalpy for internal flow too.
- Implements slip boundary condition on the north boundary of the internal
flow.
- Reorganizes the source code for the internal flow. Changes in the code did not changed the solution.
- Adjusts the extrapolation of p to fictitious corners in the internal flow.

ToDo
====
In order to make the code problem-independent, the following modifications are necessary:
* Mesh must be generated outside MACH-2D
* Boundary conditions must be defined in the configuration file

