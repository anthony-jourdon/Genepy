.. 
  ====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  Genepy
  filename: usage.rst

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


Usage
=====

.. _installation:

Installation
------------
First, this module requires python >= 3.10  and the following 
dependencies that can all be installed using pip

- numpy
- sympy
- matplotlib
- pyvista

To be able to import genepy module its location needs to be known by 
the environment variable **$PYTHONPATH**. 
To do so, move yourself in Genepy and type

.. code-block:: console

  source genepy.sh

This needs to be done every time a new terminal or python session is opened.
Another possibility is to add it to your bashrc or zshrc but it may be erased by python virtual environment when loaded.

Once your **$PYTHONPATH** has been appended check that it works correctly by typing

.. code-block:: python

  import genepy

in a python terminal. If no error is raised, the module is correctly installed.

Building a model
----------------

.. toctree:: 
  :maxdepth: 3

  gmsh_tuto

.. warning:: 
  In each examples, the following code blocks may re-use variables defined in previous blocks.
  To ensure that the code runs correctly, 
  it is recommended to run the code blocks in the order they are presented.

.. toctree:: 
  :maxdepth: 3

  examples/simple_model
  examples/oblique_model
  examples/strikeslip_model
  examples/time_dependant_model
  examples/ale_model

Run pTatin3d
------------

.. warning:: This section **does not** cover the installation of `pTatin3d`_.

  All commands are given to run in serial (1 MPI rank) and using 
  the standard bash command line arguments. For a parallel run on a HPC machine you need 
  to refer to the machine's documentation.

  In the following examples, 
  the environment variable :code:`PETSC_ARCH` is assumed to be known.

To run the model, you need to have `pTatin3d`_ installed on your system.
Once your model is built and the options file is saved, 
you can run the model using the commands presented below.

Compute initial topography
.......................... 
If your problem involves a density distribution that should produce non-zero topography,
`pTatin3d`_ provides an automatic method to compute an initial isostatic topography.
Options related to this problem are provided by default and can be changed using the 
corresponding keywords arguments found in the :py:class:`Model <genepy.Model>` class.

The option

.. code-block:: bash

  -model_GENE3D_isostatic_density_ref 3300
  
indicate the reference density in kg.m\ :sup:`-3` considered to compute the isostatic equilibrium and the option

.. code-block:: bash

  -model_GENE3D_isostatic_depth -40e3

indicate the depth at which the compensation should be computed.

.. note:: 
  As a rule of thumb from experiments, the compensation depth should be chosen
  near the transition from lower densities to the reference density.
  As an example for lithospheric models, the approximate Moho depth is a decent candidate.  

To compute this topography run the following:

.. code-block:: bash

  $PETSC_ARCH/bin/test_ptatin_driver_pressure_poisson.app -options_file path_to_file.opts -run -isostatic_remesh

It will write a file named ``isostatic_displacement.pbvec`` that will be used by the next driver
ran to adjust the topography, therefore to verify the generated topography you need to run another 
driver among the ones presented below.

Running initial conditions driver
.................................
Before running an actual solve, it is good practice to first run the initial conditions 
driver of `pTatin3d`_ to verify that the Stokes boundary conditions, the initial geometry
and the potential initial plastic strain are correctly defined.

.. note::
  If the viscosity type requested is non-linear and depends on the velocity and temperature, 
  the viscosity may not be correct because the velocity, pressure and temperature fields 
  have not been solved for.

.. code-block:: bash

  $PETSC_ARCH/bin/ptatin_driver_ic.app -options_file path_to_file.opts

By default, the following options are added to the options file.

.. code-block:: bash

  -model_GENE3D_output_markers
  -model_GENE3D_bc_debug

Before running a large job you should remove them from your options file to avoid 
the flood of the filesystem and standard output.

Running linear driver
.....................
If the problem is linear i.e., the viscosities are viscous linear you can run 

.. code-block:: bash

  $PETSC_ARCH/bin/ptatin_driver_linear_ts.app -options_file path_to_file.opts

Computing steady-state temperature
...................................
If your problem involves temperature, you can compute the initial temperature distribution
using a steady-state solution of the heat equation.
`pTatin3d`_ provides a driver to compute this solution.
Run:

.. code-block:: bash

  $PETSC_ARCH/bin/test_steady_state_diffusion_solve_TFV.app -options_file path_to_file.opts


.. warning:: 
  If your problem involves the asthenosphere, to produce a constant vertical temperature variation in the asthenosphere

  .. math:: 
    \frac{\partial T}{\partial y} = c

  i.e., a linear temperature distribution in the asthenosphere, you need to provide a high conductivity value to your asthensophere.
  However, be careful to set back a reasonable value for the conductivity before running the time dependant problem. 

This will write a file named ``temperature_steady.pbvec`` and if your options file 
contains the option (default):

.. code-block:: bash

  -view_ic

it will also output a file named ``T_steady.vts`` that contains the solution. 

Running non-linear driver with checkpointing
............................................
Finally, after computing the initial topography (if required) and 
initial temperature distribution, 
to run a non-linear problem with checkpointing capabilities you can run 

.. code-block:: bash

  $PETSC_ARCH/bin/test_ptatin_driver_checkpoint_fv.app -options_file path_to_file.opts -init
  $PETSC_ARCH/bin/test_ptatin_driver_checkpoint_fv.app -options_file path_to_file.opts -run
