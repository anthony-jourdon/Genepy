.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  Genepy
  filename: boundary_conditions.rst

  This file is part of Genepy.

  Genepy is free software: you can redistribute it and/or modify it under the terms 
  of the GNU General Public License as published by the Free Software Foundation, either 
  version 3 of the License, or any later version.

  Genepy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with Genepy. 
  If not, see <https://www.gnu.org/licenses/>.
  ====================================================================================================

.. _pTatin3d: https://github.com/laetitialp/ptatin-gene
.. _PETSc: https://petsc.org

Boundary conditions
===================
Contains a class to evaluate symbolic and numeric linear velocity functions and their derivatives 
and classes to generate options for input files of `pTatin3d <https://github.com/laetitialp/ptatin-gene>`_.   

Velocity
--------
This module contains the class evaluating symbolic and numeric velocity field.

.. autoclass:: genepy.Velocity
    :members:

Linear velocity function
........................

.. autoclass:: genepy.VelocityLinear
    :members:

Time dependant velocity function
................................

.. autoclass:: genepy.VelocityTimeDependant
    :members:

.. autoclass:: genepy.VelocityInversion
    :members:

Options generation
------------------

Bounary conditions wrapper
..........................

This class is a wrapper to generate options of the boundary conditions 
for input files of `pTatin3d`_.

.. autoclass:: genepy.ModelBCs
    :members:


Stokes boundary conditions
.......................... 
This class is parent of all Stokes boundary conditions classes:

  - :py:class:`Dirichlet <genepy.boundary_conditions.dirichlet.Dirichlet>`
  - :py:class:`NavierSlip <genepy.boundary_conditions.navierslip.NavierSlip>`
  - :py:class:`Neumann <genepy.boundary_conditions.neumann.Neumann>`

.. autoclass:: genepy.StokesBoundaryCondition
    :members:

Dirichlet
.......... 
Classes to generate options for Dirichlet boundary conditions.
Inherits from class :py:class:`StokesBoundaryCondition <genepy.boundary_conditions.bcs.StokesBoundaryCondition>`.

.. autoclass:: genepy.Dirichlet
    :members:

Navier-slip
...........
Class to generate options for Navier-slip boundary conditions.
Inherits from class :py:class:`StokesBoundaryCondition <genepy.boundary_conditions.bcs.StokesBoundaryCondition>`.

.. autoclass:: genepy.NavierSlip
    :members:

Neumann
........
Class to generate options for Neumann boundary conditions.
Inherits from class :py:class:`StokesBoundaryCondition <genepy.boundary_conditions.bcs.StokesBoundaryCondition>`.

.. autoclass:: genepy.Neumann
    :members:
