.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  bc-pre-processing
  filename: initial_conditions.rst

  This file is part of bc-pre-processing.

  bc-pre-processing is free software: you can redistribute it and/or modify it under the terms 
  of the GNU General Public License as published by the Free Software Foundation, either 
  version 3 of the License, or any later version.

  bc-pre-processing is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with bc-pre-processing. 
  If not, see <https://www.gnu.org/licenses/>.
  ====================================================================================================

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

.. autoclass:: genepy.Domain
    :members:


Mesh refinement
---------------
This module contains the class describing the mesh refinement.

.. autoclass:: genepy.MeshRefinement
    :members:

Rotation
--------
This module contains the class to perform rotations of single vectors, 
vector fields and referential in 2D and 3D.

.. autoclass:: genepy.Rotation
    :members:

Gaussian
--------
This module contains the class to evaluate gaussian distributions of a field in 2D.
It is generally used to define the initial strain distribution to place weak zones in the domain.

.. autoclass:: genepy.Gaussian
    :members:

ICs pTatin3d options generation
-------------------------------
This module contains the class to generate the options 
for the initial conditions of a 3D model running GENE3D 
in `pTatin3d`_.

.. autoclass:: genepy.InitialConditions
    :members:
