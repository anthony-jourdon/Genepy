# Boundary conditions pre-processing
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