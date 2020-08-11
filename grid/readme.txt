Folder "grid" contains modules of data and procedures related to the grid.

-       mod_grid_data.f90: contains DATA and type definitions related to the grid.
- mod_grid_procedures.f90: contains PROCEDURES to operate grid data.
-            mod_grid.f90: creates an interface to the main program to calculate the grid and related
                           metrics. Uses PROCEDURES from mod_grid_procedures to operate DATA from mod_grid_data.

