# Boundary conditions pre-processing
## Presentation
`bcpy` is a python module designed to evaluate symbolic (mathematic) expression to build analytical velocity functions varying in space and time, define the rheological parameters of a long-term geodynamic model, and generate input file for [pTatin3d](https://github.com/laetitialp/ptatin-gene).
It is a pre-processing tool.

This module can:
- evaluate and print mathematical expression for the velocity and initial plastic strain distribution 
- attribute rheological parameters to regions identified by a tag (integer value)
- handle simple mesh refinement for a structured mesh using linear interpolation
- generate options file for pTatin3d simulations

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

## Documentation
To build the documentation of this package go to the root directory of the repository and type

```
sphinx-build -M html docs/source/ docs/build/
```

in docs/build/ you should find index.hmtl, open it with your web browser and navigate through the documentation!