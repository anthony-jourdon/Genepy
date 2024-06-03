.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  bc-pre-processing
  filename: material_parameters.rst

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

Material parameters
===================

Wrapper for material parameters
-------------------------------

Region
......
This class is a wrapper for material parameters. 
It is used to set all material parameters of a given region at once.

.. autoclass:: bcpy.Region
    :members:

Model regions
.............

.. autoclass:: bcpy.ModelRegions
    :members:

Parent class for all material parameters.
-----------------------------------------

.. autoclass:: bcpy.MaterialConstants
    :members:

.. _density:

Density
-------

Constant density
................ 

.. autoclass:: bcpy.DensityConstant
    :members:

Boussinesq approximation
........................ 

.. autoclass:: bcpy.DensityBoussinesq
    :members:

Thermodynamic table
...................

.. autoclass:: bcpy.DensityTable
    :members:

.. _viscosity:

Viscosity
---------

Constant viscosity
..................

.. autoclass:: bcpy.ViscosityConstant
    :members:

Frank-Kamenetskii viscosity
............................

.. autoclass:: bcpy.ViscosityFrankK
    :members:

Depth-dependent viscosity
..........................

.. autoclass:: bcpy.ViscosityZ
    :members:

Arrhenius viscosity
....................

.. autoclass:: bcpy.ViscosityArrhenius
    :members:

.. autoclass:: bcpy.ViscosityArrhenius2
    :members:

.. autoclass:: bcpy.ViscosityArrheniusDislDiff
    :members:

.. _plasticity:

Pasticity
---------

No plasticity
.............

.. autoclass:: bcpy.PlasticNone
    :members:

VonMises plasticity
....................

.. autoclass:: bcpy.PlasticMises
    :members:

Drucker-Prager plasticity
..........................

.. autoclass:: bcpy.PlasticDruckerPrager
    :members:

.. _softening:

Plastic strain softening
------------------------

No softening
.............

.. autoclass:: bcpy.SofteningNone
    :members:

Linear softening
................

.. autoclass:: bcpy.SofteningLinear
    :members:

Exponential softening
.....................

.. autoclass:: bcpy.SofteningExponential
    :members:

.. _energy:

Energy
------

.. autoclass:: bcpy.Energy
    :members:
