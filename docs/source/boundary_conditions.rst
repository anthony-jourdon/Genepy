.. _pTatin3d: https://github.com/laetitialp/ptatin-gene
.. _PETSc: https://petsc.org

Boundary conditions
===================
Contains a class to evaluate symbolic and numeric linear velocity functions and their derivatives 
and classes to generate options for input files of `pTatin3d <https://github.com/laetitialp/ptatin-gene>`_.   

Velocity
--------
This module contains the class evaluating symbolic and numeric velocity field.

.. autoclass:: bcpy.Velocity
    :members:

Options generation
------------------

Bounary conditions wrapper
..........................

This class is a wrapper to generate options of the boundary conditions 
for input files of `pTatin3d`_.

.. autoclass:: bcpy.ModelBCs
    :members:


Stokes boundary conditions
.......................... 
This class is parent of all Stokes boundary conditions classes:

  - :py:class:`Dirichlet <bcpy.boundary_conditions.dirichlet.Dirichlet>`
  - :py:class:`NavierSlip <bcpy.boundary_conditions.navierslip.NavierSlip>`
  - :py:class:`Neumann <bcpy.boundary_conditions.neumann.Neumann>`

.. autoclass:: bcpy.StokesBoundaryCondition
    :members:

Dirichlet
.......... 
Classes to generate options for Dirichlet boundary conditions.
Inherits from class :py:class:`StokesBoundaryCondition <bcpy.boundary_conditions.bcs.StokesBoundaryCondition>`.

.. autoclass:: bcpy.Dirichlet
    :members:

Navier-slip
...........
Class to generate options for Navier-slip boundary conditions.
Inherits from class :py:class:`StokesBoundaryCondition <bcpy.boundary_conditions.bcs.StokesBoundaryCondition>`.

.. autoclass:: bcpy.NavierSlip
    :members:

Neumann
........
Class to generate options for Neumann boundary conditions.
Inherits from class :py:class:`StokesBoundaryCondition <bcpy.boundary_conditions.bcs.StokesBoundaryCondition>`.

.. autoclass:: bcpy.Neumann
    :members: