# Boundary conditions pre-processing
## Presentation
`bcpy` is a python module designed to evaluate symbolic (mathematic) expression to build analytical velocity functions varying in space and time.

Its first goal is to be used as a pre-processing tool to get boundary conditions for long-term geodynamic thermo-mechanical models, in particular [pTatin3d](https://bitbucket.org/dmay/ptatin-total-dev).

This module can:
- evaluate and print mathematical expression for the velocity 
- allow visualisation of the velocity field using matplotlib if steady-state or paraview (.vts + timeseries.pvd files) for time varying velocities

Some examples can be found in `scripts`  subdirectory.

## Installation
To be able to import bcpy module its location needs to be known by `PYTHONPATH`.
To do so, move yourself in bc-pre-processing and type

`source bc-pre-proc.sh`

This needs to be done every time a new terminal or python session is opened.
Another possibility is to add it to your bashrc or zshrc but it may be erased by python virtual environment when loaded.

Once your `PYTHONPATH` has been appended check that it works correctly with

```
$ python
>> import bcpy
```

This module requires the following dependencies that can all be installed using pip
- numpy
- sympy
- matplotlib
- pyvista