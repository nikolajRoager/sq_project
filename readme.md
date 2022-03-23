Sketch of particles in electromagnetic fields
--------
This project simulates one or more classical particles moving in electric and magnetic fields.

by default, the makefile builds the project and runs the project with the setup file `assets\setup.txt` and saves terminal log to `log/log_default.txt` (mainly for debugging) and the output to the folder `out\out_default`. output and log can be cleared with `make clean_log` or `make clean_out`. (Saving new outputs, labelled by the time and data, every time make was run was tried but removed, as it created a lot of excess data)

To display the simulation result, use the provided python file `python display_data.py out_folder (args)`, `out_folder` must be the folder where the simulation put the output, args is how to display, the options are `--3D` (default), which shows the fields and particle in a 3D plot. 3D vector plots are, however, not very easy to read.

A more readable alternative is to use the arugment `--x number` (or `--y` or `--z`), this only shows a slice of the field at `x,y or z = whatever` (In whatever unit was defined in the setup file). This display option flattens the path taken by the particle. 2D displays essentially throw away whatever dimension is deemed less important.

The setup file defines what field is used, how many particles are used (and where they start) how long the simulation runs and what resolution is used, and what to output for display.
