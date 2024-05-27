.. _pTatin3d: https://github.com/laetitialp/ptatin-gene
.. _PETSc: https://petsc.org

Initial conditions
==================
Contains the classes to define the initial conditions such as the 
geometry of the physical domain, the mesh refinement or an initial field with 
a gaussian distribution.
It generates symbolic expressions of the coordinates and initial field distribution.

Domain
------
In `pTatin3d`_, the physical domain where the simulation is performed
is defined by the coordinate system showed in the figure below.

.. image:: figures/coordinate_system.png
    :width: 400
    :align: center


The following code describes the physical domain.
While most of the usage of this class is for 3 dimensional domains,
it can also be used for 2 dimensional and 1 dimensional domains.

.. autoclass:: bcpy.Domain
    :members:


Mesh refinement
---------------
This module contains the class describing the mesh refinement.

.. autoclass:: bcpy.MeshRefinement
    :members:

Rotation
--------
This module contains the class to perform rotations of single vectors, 
vector fields and referential in 2D and 3D.

.. autoclass:: bcpy.Rotation
    :members:

Gaussian
--------
This module contains the class to evaluate gaussian distributions of a field in 2D.
It is generally used to define the initial strain distribution to place weak zones in the domain.

.. autoclass:: bcpy.Gaussian
    :members:

ICs pTatin3d options generation
-------------------------------
This module contains the class to generate the options 
for the initial conditions of a 3D model running GENE3D 
in `pTatin3d`_.

.. autoclass:: bcpy.InitialConditions
    :members:
