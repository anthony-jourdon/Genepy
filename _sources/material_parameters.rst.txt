.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  Genepy
  filename: material_parameters.rst

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

Material parameters
===================

Wrapper for material parameters
-------------------------------

Region
......
This class is a wrapper for material parameters. 
It is used to set all material parameters of a given region at once.

.. autoclass:: genepy.Region
    :members:

Model regions
.............

.. autoclass:: genepy.ModelRegions
    :members:

Parent class for all material parameters.
-----------------------------------------

.. autoclass:: genepy.MaterialConstants
    :members:

.. _density:

Density
-------

Constant density
................ 

.. autoclass:: genepy.DensityConstant
    :members:

Boussinesq approximation
........................ 

.. autoclass:: genepy.DensityBoussinesq
    :members:

Thermodynamic table
...................

.. autoclass:: genepy.DensityTable
    :members:

.. _viscosity:

Viscosity
---------

Constant viscosity
..................

.. autoclass:: genepy.ViscosityConstant
    :members:

Frank-Kamenetskii viscosity
............................

.. autoclass:: genepy.ViscosityFrankK
    :members:

Depth-dependent viscosity
..........................

.. autoclass:: genepy.ViscosityZ
    :members:

Arrhenius viscosity
....................

.. autoclass:: genepy.ViscosityArrhenius
    :members:

.. autoclass:: genepy.ViscosityArrhenius2
    :members:

.. autoclass:: genepy.ViscosityArrheniusDislDiff
    :members:

.. _plasticity:

Pasticity
---------

No plasticity
.............

.. autoclass:: genepy.PlasticNone
    :members:

VonMises plasticity
....................

.. autoclass:: genepy.PlasticMises
    :members:

Drucker-Prager plasticity
..........................

.. autoclass:: genepy.PlasticDruckerPrager
    :members:

.. _softening:

Plastic strain softening
------------------------

No softening
.............

.. autoclass:: genepy.SofteningNone
    :members:

Linear softening
................

.. autoclass:: genepy.SofteningLinear
    :members:

Exponential softening
.....................

.. autoclass:: genepy.SofteningExponential
    :members:

.. _energy:

Energy
------

.. autoclass:: genepy.Energy
    :members:
