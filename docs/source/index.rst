.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  Genepy
  filename: index.rst

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

.. Genepy documentation master file, created by
   sphinx-quickstart on Thu May 23 17:19:03 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Genepy's documentation!
==================================
`genepy` is a python module designed to evaluate symbolic (mathematic) 
expression to build analytical velocity functions varying in space and time, 
define the rheological parameters of a long-term geodynamic model, 
and generate input file for `pTatin3d`_.
It is a pre-processing tool.

This module can:

- evaluate and print mathematical expression for the velocity and initial plastic strain distribution 

- attribute rheological parameters to regions identified by a tag (integer value)

- handle simple mesh refinement for a structured mesh using linear interpolation

- generate options file for pTatin3d simulations

Some examples can be found in *scripts*  subdirectory.
Check out the :doc:`usage` section for further information, including how to
:ref:`install <installation>` the project.

.. note::
    This module is still under development.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Contents
========

.. toctree:: 

   usage
   genepy_docs
