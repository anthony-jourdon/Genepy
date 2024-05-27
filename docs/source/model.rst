.. _pTatin3d: https://github.com/laetitialp/ptatin-gene
.. _PETSc: https://petsc.org

Model
=====
This class generates the options for `pTatin3d`_ model using
the format of the `PETSc`_ library.
It requires to have generated 

  - the physical :py:class:`domain <bcpy.Domain>`
  - the :py:class:`velocity <bcpy.Velocity>` field function.
  - the :py:class:`initial conditions <bcpy.InitialConditions>`.
  - the :py:class:`boundary conditions <bcpy.ModelBCs>`.
  - the :py:class:`material parameters <bcpy.ModelRegions>`.

.. autoclass:: bcpy.Model
    :members:

Markers management
==================

.. autoclass:: bcpy.MarkersManagement
    :members:

Passive tracers
===============

Pswarm
------

.. autoclass:: bcpy.Pswarm
    :members:

Fill entire domain
------------------

.. autoclass:: bcpy.PswarmFillDomain
    :members:

Fill domain within a bounding box
---------------------------------

.. autoclass:: bcpy.PswarmFillDomainWithinBoundingBox
    :members:

Fill domain within a bounding box with a given number of particles
------------------------------------------------------------------

.. autoclass:: bcpy.PswarmFillBox
    :members:

Place only selected tracers
---------------------------

.. autoclass:: bcpy.PswarmFromUserList
    :members:

Surface processes
=================

.. autoclass:: bcpy.SPMDiffusion
    :members: