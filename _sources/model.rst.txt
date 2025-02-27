.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  Genepy
  filename: model.rst

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


Model
=====
This class generates the options for `pTatin3d`_ model using
the format of the `PETSc`_ library.
It requires to have generated 

  - the physical :py:class:`domain <genepy.Domain>`
  - the :py:class:`velocity <genepy.VelocityLinear>` field function.
  - the :py:class:`initial conditions <genepy.InitialConditions>`.
  - the :py:class:`boundary conditions <genepy.ModelBCs>`.
  - the :py:class:`material parameters <genepy.ModelRegions>`.

.. autoclass:: genepy.Model
    :members:

Markers management
==================

.. autoclass:: genepy.MarkersManagement
    :members:

Passive tracers
===============

Pswarm
------

.. autoclass:: genepy.Pswarm
    :members:

Fill entire domain
------------------

.. autoclass:: genepy.PswarmFillDomain
    :members:

Fill domain within a bounding box
---------------------------------

.. autoclass:: genepy.PswarmFillDomainWithinBoundingBox
    :members:

Fill domain within a bounding box with a given number of particles
------------------------------------------------------------------

.. autoclass:: genepy.PswarmFillBox
    :members:

Place only selected tracers
---------------------------

.. autoclass:: genepy.PswarmFromUserList
    :members:

Surface processes
=================

.. autoclass:: genepy.SPMDiffusion
    :members:

.. autoclass:: genepy.SPMDiffusionBaselvl
    :members: