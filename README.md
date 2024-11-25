# Genepy
## Presentation
`genepy` is a python module designed to evaluate symbolic (mathematic) expression to build analytical velocity functions varying in space and time, define the rheological parameters of a long-term geodynamic model, and generate input file for [pTatin3d](https://github.com/laetitialp/ptatin-gene).
It is a pre-processing tool.

This module can:
- evaluate and print mathematical expression for the velocity and initial plastic strain distribution 
- attribute rheological parameters to regions identified by a tag (integer value)
- handle simple mesh refinement for a structured mesh using linear interpolation
- generate options file for pTatin3d simulations

Some examples can be found in `scripts`  subdirectory.

## Installation
This module requires python >= 3.10 and the following dependencies that can all be installed using pip
- numpy
- sympy
- matplotlib
- pyvista
- gmsh

To be able to import genepy module its location needs to be known by `PYTHONPATH`.
To do so, move yourself in Genepy and type

`source genepy.sh`

This needs to be done every time a new terminal or python session is opened.
Another possibility is to add it to your bashrc or zshrc but it may be erased by python virtual environment when loaded.

Once your `PYTHONPATH` has been appended check that it works correctly with

```
$ python
>> import genepy
```

## Documentation
Access documentation with this [link](https://anthony-jourdon.github.io/Genepy/).

To update the documentation install `ghp-import` python package found [here](https://github.com/c-w/ghp-import) and use the following command

```
ghp-import docs/build/html -o -p -n
```
