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

.. _pTatin3d: https://github.com/laetitialp/ptatin-gene

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

  source bc-pre-proc.sh

This needs to be done every time a new terminal or python session is opened.
Another possibility is to add it to your bashrc or zshrc but it may be erased by python virtual environment when loaded.

Once your **$PYTHONPATH** has been appended check that it works correctly by typing

.. code-block:: python

  import genepy

in a python terminal. If no error is raised, the module is correctly installed.

Building a model
----------------
.. warning:: 
  In each examples, the following code blocks may re-use variables defined in previous blocks.
  To ensure that the code runs correctly, 
  it is recommended to run the code blocks in the order they are presented.

Example: simple model, linear viscous rheology
..............................................
This example shows how to build a very simple model using default linear viscous rheology 
and standard Dirichlet type boundary conditions for a velocity field imposing shortening in
:math:`z` direction as shown by the figure below.
Because the viscosity chosen in this example does not depend on temperature,
the thermal part of the model is not included.

.. image:: figures/compression.png
   :width: 200
   :align: center


1. Create a domain
~~~~~~~~~~~~~~~~~~~
We define a 3D :py:class:`domain <genepy.Domain>` :math:`\Omega = [0,600]\times[-250,0]\times[0,300]` km\ :sup:`3`
:math:`\in \mathbb R^3` discretized by a regular grid of 9x9x9 nodes.

.. code-block:: python

  import os
  import numpy as np
  import genepy as bp

  # 3D domain
  dimensions = 3
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = bp.Domain(dimensions,O,L,n)

2. Velocity function
~~~~~~~~~~~~~~~~~~~~
We define a simple orthogonal shortening :py:class:`velocity <genepy.Velocity>` in the :math:`z` direction.

.. code-block:: python

  # velocity
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "compression"                  # extension or compression
  # Create Velocity class instance
  BCs = bp.Velocity(Domain,u_norm,u_dir,u_type)

  # Evaluate the velocity and its derivatives
  u,grad_u = BCs.evaluate_velocity_and_gradient_symbolic() # symbolic

.. note:: In this example, the derivatives of the velocity are not used.

3. Initial conditions
~~~~~~~~~~~~~~~~~~~~~
In this example we do not impose any initial plastic strain value nor mesh refinement.
Therefore the :py:class:`initial conditions <genepy.InitialConditions>` are only the Domain and the velocity function.
They will be used to generate the options for `pTatin3d`_ model.

.. code-block:: python

  # Initial conditions
  model_ics = bp.InitialConditions(Domain,u)

4. Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~
Because the imposed velocity is orthogonal to the boundary we can define the 
velocity boundary conditions using :py:class:`Dirichlet <genepy.Dirichlet>` type 
:py:class:`boundary conditions <genepy.ModelBCs>`.

.. note:: 
  In the following example a path to the mesh files describing the boundaries is provided.
  These mesh files are located in ``"ptatin-gene/src/models/gene3d/examples"``.
  You can modify the ``root`` variable to match the location of the mesh files on your system 
  or remove that part of the code if you do not have access to these files.
  Note however that `pTatin3d`_ requires mesh files to define the boundaries.

Details on the methods used to define the boundary conditions can be found in the
:doc:`boundary conditions <boundary_conditions>` section.

.. code-block:: python

  # boundary conditions
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  u_bcs = [
    bp.Dirichlet(23,"Zmax",["z"],u, mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")), # orthogonal shortening
    bp.Dirichlet(37,"Zmin",["z"],u, mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")), # orthogonal shortening
    bp.Dirichlet(32,"Xmax",["x"],u, mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")), # free-slip
    bp.Dirichlet(14,"Xmin",["x"],u, mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")), # free-slip
    bp.DirichletUdotN(33,"Bottom",  mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")), # basal outflow
  ]
  # collect all boundary conditions
  model_bcs = bp.ModelBCs(u_bcs)

5. Material parameters
~~~~~~~~~~~~~~~~~~~~~~
Next we define the material properties of each :py:class:`Region <genepy.Region>` and 
gather them all in a :py:class:`ModelRegions <genepy.ModelRegions>` class instance.
In this example we use the default values for **all regions**:

- :py:class:`Constant viscosity <genepy.ViscosityConstant>` of :math:`10^{22}` Pa.s.
- :py:class:`Constant density <genepy.DensityConstant>` of :math:`3300` kg.m\ :sup:`-3`.
- :py:class:`No plasticity <genepy.PlasticNone>`.
- :py:class:`No softening <genepy.SofteningNone>`.

.. code-block:: python

  regions = [
    # Upper crust
    bp.Region(38),
    # Lower crust
    bp.Region(39),
    # Lithosphere mantle
    bp.Region(40),
    # Asthenosphere
    bp.Region(41)
  ]
  model_regions = bp.ModelRegions(regions,
                                  mesh_file=os.path.join(root,"box_ptatin_md.bin"),
                                  region_file=os.path.join(root,"box_ptatin_region_cell.bin"))

6. Create the model and generate options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, we create the :py:class:`model <genepy.Model>` by gathering all the information defined previously and we save
the options to a file named ``simple_shortening_model.opts``.

.. code-block:: python

  # create class instance
  model = bp.Model(model_ics,model_regions,model_bcs)
  # write the options for ptatin3d
  with open("simple_shortening_model.opts","w") as f:
    f.write(model.options)

Example: oblique model, non-linear rheology
...........................................
In this example we build a model with an oblique velocity field to impose 
extension at 30 degrees (counter-clockwise) with respect to the :math:`z` axis 
(can be seen as north-south direction).
We use :py:class:`non-linear viscous <genepy.ViscosityArrhenius2>` rheology, 
:py:class:`Drucker-Prager plasticity <genepy.PlasticDruckerPrager>` and
a combination of :py:class:`Dirichlet <genepy.Dirichlet>` and 
:py:class:`Navier-slip <genepy.NavierSlip>` type boundary conditions.

.. image:: figures/Oblique_extension.PNG
   :width: 400
   :align: center

1. Create a domain
~~~~~~~~~~~~~~~~~~~
We define a 3D domain :math:`\Omega = [0,600]\times[-250,0]\times[0,300]` km\ :sup:`3`
:math:`\in \mathbb R^3` discretized by a regular grid of 9x9x9 nodes.

.. code-block:: python

  import os
  import numpy as np
  import genepy as bp

  # 3D domain
  dimensions = 3
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = bp.Domain(dimensions,O,L,n)

2. Velocity function
~~~~~~~~~~~~~~~~~~~~
We define an oblique extension :py:class:`velocity <genepy.Velocity>` velocity field
forming an angle of 30 degrees counter-clockwise with respect to the :math:`z` axis.
The method 
:py:meth:`evaluate_velocity_and_gradient_symbolic() <genepy.Velocity.evaluate_velocity_and_gradient_symbolic>` 
returns the symbolic expression of the velocity field and its gradient.
The method
:py:meth:`evaluate_velocity_numeric() <genepy.Velocity.evaluate_velocity_numeric>`
returns the numeric value of the velocity field evaluated at coordinates of the nodes.
The method
:py:meth:`get_velocity_orientation() <genepy.Velocity.get_velocity_orientation>`
returns the orientation of the velocity field at the boundary.

.. code-block:: python

  # velocity
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_angle = np.deg2rad(30.0)               # velocity angle \in [-pi/2, pi/2]
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "extension"                    # extension or compression
  # Create Velocity class instance
  BCs = bp.Velocity(Domain,u_norm,u_dir,u_type,u_angle)

  # Evaluate the velocity and its derivatives
  u,grad_u = BCs.evaluate_velocity_and_gradient_symbolic() # symbolic
  u_num    = BCs.evaluate_velocity_numeric()                  # numeric
  uL       = BCs.get_velocity_orientation(horizontal=True,normalize=True)

3. Initial conditions
~~~~~~~~~~~~~~~~~~~~~
In this example we do not impose any initial plastic strain value nor mesh refinement.
Therefore the :py:class:`initial conditions <genepy.InitialConditions>` 
are only the Domain and the velocity function.
They will be used to generate the options for `pTatin3d`_ model.

.. code-block:: python

  # Initial conditions
  model_ics = bp.InitialConditions(Domain,u)

4. Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~
Because the imposed velocity is oblique to the boundary we define the
velocity boundary conditions using :py:class:`Dirichlet <genepy.Dirichlet>` and
:py:class:`Navier-slip <genepy.NavierSlip>` type :py:class:`boundary conditions <genepy.ModelBCs>`.
Note that the Dirichlet conditions takes now the 2 horizontal components to impose the obliquity. 

Moreover, we will use non-linear viscosities depending of the temperature 
so we need to provide boundary conditions for the conservation of the thermal energy.

Details on the methods used to define the boundary conditions can be found in the
:doc:`boundary conditions <boundary_conditions>` section.

.. code-block:: python

  # boundary conditions
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  u_bcs = [
    bp.Dirichlet( 23,"Zmax",["x","z"],u, mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")),
    bp.Dirichlet( 37,"Zmin",["x","z"],u, mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")),
    bp.NavierSlip(32,"Xmax",grad_u,uL,   mesh_file=os.path.join(root,"box_ptatin_facet_32_mesh.bin")),
    bp.NavierSlip(14,"Xmin",grad_u,uL,   mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
    bp.DirichletUdotN(33,"Bottom",       mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")),
  ]
  # Temperature boundary conditions
  Tbcs = bp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  model_bcs = bp.ModelBCs(u_bcs,Tbcs)

5. Material parameters
~~~~~~~~~~~~~~~~~~~~~~
Next we define the material properties of each :py:class:`Region <genepy.Region>` and
gather them all in a :py:class:`ModelRegions <genepy.ModelRegions>` class instance.
In this example we use the following material types:

- :py:class:`Dislocation creep <genepy.ViscosityArrhenius2>`.
- :py:class:`Drucker-Prager <genepy.PlasticDruckerPrager>` plastic yield criterion.
- :py:class:`Linear softening <genepy.SofteningLinear>`.
- :py:class:`Boussinesq density <genepy.DensityBoussinesq>`.

.. code-block:: python

  regions = [
    # Upper crust
    bp.Region(38,                                          # region tag
              bp.DensityBoussinesq(2700.0,3.0e-5,1.0e-11), # density
              bp.ViscosityArrhenius2("Quartzite"),         # viscosity  (values from the database using rock name)
              bp.SofteningLinear(0.0,0.5),                 # softening
              bp.PlasticDruckerPrager(),                   # plasticity (default values, can be modified using the corresponding parameters)
              bp.Energy(1.5e-6,2.7)),                      # energy
    # Lower crust
    bp.Region(39,
              bp.DensityBoussinesq(density=2850.0,thermal_expansion=3.0e-5,compressibility=1.0e-11),
              bp.ViscosityArrhenius2("Anorthite",Vmol=38.0e-6),
              bp.SofteningLinear(strain_min=0.0,strain_max=0.5),
              bp.PlasticDruckerPrager(),
              bp.Energy(heat_source=0.5e-6,conductivity=2.85)),
    # Lithosphere mantle
    bp.Region(40,
              bp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              bp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.Energy(0.0,3.3)),
    # Asthenosphere
    bp.Region(41,
              bp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              bp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.Energy(0.0,3.3))
  ]
  model_regions = bp.ModelRegions(regions,
                                  mesh_file=os.path.join(root,"box_ptatin_md.bin"),
                                  region_file=os.path.join(root,"box_ptatin_region_cell.bin"))

6. Create the model and generate options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, we create the :py:class:`model <genepy.Model>` by gathering all the information defined previously and we save
the options to a file named ``oblique_extension_model.opts``.

.. code-block:: python

  # create class instance
  model = bp.Model(model_ics,model_regions,model_bcs)
  # write the options for ptatin3d
  with open("oblique_extension_model.opts","w") as f:
    f.write(model.options)

Example: strike-slip model, rotated velocity field and mesh refinement
......................................................................
This example will build a 3D model with vertical 
:py:class:`mesh refinement <genepy.MeshRefinement>` 
and a strike-slip velocity field 
:py:class:`rotated <genepy.Rotation>` 
by 15 degrees as showed in the figure below.
In addition, 2 :py:class:`gaussian <genepy.Gaussian>` weak zones are added to the initial conditions of the model 

.. image:: figures/Strike_slip-01.png
   :width: 600
   :align: center

1. Create a domain
~~~~~~~~~~~~~~~~~~~
We define a 3D :py:class:`Domain <genepy.Domain>` :math:`\Omega = [0,600]\times[-250,0]\times[0,300]` km\ :sup:`3` 
:math:`\in \mathbb R^3` discretized by a regular grid of 9x9x9 nodes. 

.. code-block:: python

  import os
  import numpy as np
  import genepy as bp

  # 3D domain
  dimensions = 3
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = bp.Domain(dimensions,O,L,n)

2. Mesh refinement
~~~~~~~~~~~~~~~~~~
In this step we :py:class:`refine the mesh <genepy.MeshRefinement>` 
in the vertical direction (:math:`y`) using linear interpolation.
Note however that the mesh refinement can be done in any direction following the same pattern.

.. code-block:: python

  # Define refinement parameters in a dictionary
  refinement = {"y": # direction of refinement
                    {"x_initial": np.array([-250,-180,-87.5,0], dtype=np.float64)*1e3, # xp
                     "x_refined": np.array([-250,-50,-16.25,0], dtype=np.float64)*1e3} # f(xp)
               }
  # Create MeshRefinement class instance
  MshRef = bp.MeshRefinement(Domain,refinement)
  # Refine the mesh
  MshRef.refine()

3. Rotation
~~~~~~~~~~~
To rotate the velocity field we first need to 
set the parameters of this :py:class:`rotation <genepy.Rotation>`.
In this example we perform a rotation of 15 degrees 
clockwise around the :math:`y` axis.

.. code-block:: python

  # Rotation of the referential
  r_angle = np.deg2rad(-15.0)                   # Rotation angle \in [-pi, pi]
  axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
  # Create instance of Rotation class
  Rotation = bp.Rotation(dimensions,r_angle,axis)

4. Velocity field
~~~~~~~~~~~~~~~~~
Next, we create a strike-slip velocity field with a norm of 1 cm.a\ :sup:`-1`.
The method 
:py:meth:`evaluate_velocity_and_gradient_symbolic() <genepy.Velocity.evaluate_velocity_and_gradient_symbolic>` 
returns the symbolic expression of the velocity field and its gradient.
The method
:py:meth:`evaluate_velocity_numeric() <genepy.Velocity.evaluate_velocity_numeric>`
returns the numeric value of the velocity field evaluated at coordinates of the nodes.
The method
:py:meth:`get_velocity_orientation() <genepy.Velocity.get_velocity_orientation>`
returns the orientation of the velocity field at the boundary.

.. note:: The rotation of the velocity field is handled inside the velocity function evaluation
  and does not require any additional step.

.. code-block:: python

  # velocity function parameters
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_angle = np.deg2rad(90.0)               # velocity angle \in [-pi/2, pi/2]
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "extension"                    # extension or compression, defines the sign
  # Create velocity class instance
  BCs = bp.Velocity(Domain,u_norm,u_dir,u_type,u_angle,Rotation)

  # Evaluate the velocity function and its derivatives
  u,grad_u = BCs.evaluate_velocity_and_gradient_symbolic() # symbolic
  u_num    = BCs.evaluate_velocity_numeric()                  # numeric
  # Get the orientation of the vectors at boundary (horizontal removes the vertical component)
  uL       = BCs.get_velocity_orientation(horizontal=True,normalize=True)

5. Define gaussian weak zones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this exemple we define two :py:class:`gaussian <genepy.Gaussian>` weak zones.
We provide the parameters of the gaussians and their position in the domain.

.. note:: 
  In this example we rotate the velocity field by 15 degrees.
  Therefore we also rotate the gaussians by 15 degrees.
  This is achieved by passing the 
  :py:class:`Rotation <genepy.Rotation>` class instance to the 
  :py:class:`Gaussian <genepy.Gaussian>` class constructor.

.. code-block:: python

  # gaussian weak zones
  ng = np.int32(2) # number of gaussians
  A  = np.array([1.0, 1.0],dtype=np.float64) # amplitude (will be multiplied by a random number between 0 and 1 to generate noise in the model)
  # coefficients for the shape of the gaussians
  coeff = 0.5 * 6.0e-5**2
  a = np.array([coeff, coeff], dtype=np.float64)
  b = np.array([0.0, 0.0],     dtype=np.float64)
  c = np.array([coeff, coeff], dtype=np.float64)
  # position of the centre of the gaussians
  dz    = 25.0e3                            # distance from the domain centre in z direction
  angle = np.deg2rad(83.0)                  # angle between the x-axis and the line that passes through the centre of the domain and the centre of the gaussian
  domain_centre = 0.5*(Domain.O + Domain.L) # centre of the domain
  
  x0 = np.zeros(shape=(ng), dtype=np.float64)
  # centre of the gaussian in z direction
  z0 = np.array([domain_centre[2] - dz, 
                 domain_centre[2] + dz], dtype=np.float64) 
  # centre of the gaussian in x direction
  x0[0] = bp.utils.x_centre_from_angle(z0[0],angle,(domain_centre[0],domain_centre[2])) 
  x0[1] = bp.utils.x_centre_from_angle(z0[1],angle,(domain_centre[0],domain_centre[2]))
  # Create instance of Gaussian class
  Gaussian = bp.Gaussian(MshRef,Rotation,ng,A,a,b,c,x0,z0)
  # Evaluate symbolic expression and numerical values of the gaussians
  Gaussian.evaluate_gaussians()

6. Initial conditions
~~~~~~~~~~~~~~~~~~~~~
Gather the information defined previously to generate the options for the initial conditions.

.. code-block:: python

  # Initial conditions
  model_ics = bp.InitialConditions(Domain,u,mesh_refinement=MshRef,initial_strain=IniStrain)

7. Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~
Gather the velocity field information and indicate the type of boundary conditions required
to generate the options for the boundary conditions.

Details on the methods used to define the boundary conditions can be found in the
:doc:`boundary conditions <boundary_conditions>` section.

.. code-block:: python

  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  u_bcs = [
    bp.Dirichlet(tag=23,name="Zmax",components=["x","z"],velocity=u,mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")),
    bp.Dirichlet(37,"Zmin",["x","z"],u,mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")),
    bp.NavierSlip(tag=32,name="Xmax",grad_u=grad_u,u_orientation=uL,mesh_file=os.path.join(root,"box_ptatin_facet_32_mesh.bin")),
    bp.NavierSlip(14,"Xmin",grad_u,uL,mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
    bp.DirichletUdotN(33,"Bottom",mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")),
  ]
  # Temperature boundary conditions
  Tbcs = bp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  model_bcs = bp.ModelBCs(u_bcs,Tbcs)

8. Material parameters
~~~~~~~~~~~~~~~~~~~~~~
Next we define the material properties (mechanical and thermal) of the different
regions of the model.
For each region, a set of parameters is defined using the corresponding classes.
The details on the methods can be found in the
:doc:`material parameters <material_parameters>` section.

.. code-block:: python

  # Define the material parameters for the model as a list of Region objects
  regions = [
    # Upper crust
    bp.Region(38,                                          # region tag
              bp.DensityBoussinesq(2700.0,3.0e-5,1.0e-11), # density
              bp.ViscosityArrhenius2("Quartzite"),         # viscosity  (values from the database using rock name)
              bp.SofteningLinear(0.0,0.5),                 # softening
              bp.PlasticDruckerPrager(),                   # plasticity (default values, can be modified using the corresponding parameters)
              bp.Energy(1.5e-6,2.7)),                      # energy
    # Lower crust
    bp.Region(39,
              bp.DensityBoussinesq(density=2850.0,thermal_expansion=3.0e-5,compressibility=1.0e-11),
              bp.ViscosityArrhenius2("Anorthite",Vmol=38.0e-6),
              bp.SofteningLinear(strain_min=0.0,strain_max=0.5),
              bp.PlasticDruckerPrager(),
              bp.Energy(heat_source=0.5e-6,conductivity=2.85)),
    # Lithosphere mantle
    bp.Region(40,
              bp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              bp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.Energy(0.0,3.3)),
    # Asthenosphere
    bp.Region(41,
              bp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              bp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.Energy(0.0,3.3))
  ]

  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  model_regions = bp.ModelRegions(regions,
                                  mesh_file=os.path.join(root,"box_ptatin_md.bin"),
                                  region_file=os.path.join(root,"box_ptatin_region_cell.bin"))

9. Add surface processes
~~~~~~~~~~~~~~~~~~~~~~~~
In this example we add :py:class:`surface processes <genepy.SPMDiffusion>`.
Surface processes are done by solving a diffusion equation. 
Here we set ``"zmin"`` and ``"zmax"`` as Dirichlet boundary conditions for the diffusion equation
and we set the diffusivity to :math:`10^6` m\ :sup:`2`.s\ :sup:`-1`.

.. code-block:: python

  # Add erosion-sedimentation with diffusion
  spm = bp.SPMDiffusion(["zmin","zmax"],diffusivity=1.0e-6)

11. Add passive tracers
~~~~~~~~~~~~~~~~~~~~~~~~
Add passive tracers to the model.
Here we define a box :math:`x \in [0, 600] \times y \in [-100, 0] \times z \in [0, 300]` km\ :sup:`3` 
of passive tracers with a layout of :math:`30 \times 5 \times 15` lagrangian markers.
We activate the tracking of the pressure and temperature fields.

.. note:: Other types of passive tracers layout can be found in the 
  :py:class:`passive tracers <genepy.Pswarm>` section.

.. code-block:: python

  # Add passive tracers
  pswarm = bp.PswarmFillBox([0.0,-100.0e3,0.0],
                            [600e3,-4.0e3,300.0e3],
                            layout=[30,5,15],
                            pressure=True,
                            temperature=True)

12.  Create the model and generate options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :py:class:`model <genepy.Model>` is created by gathering all the information defined previously.

.. code-block:: python

  # write the options for ptatin3d
  model = bp.Model(model_ics,model_regions,model_bcs,
                   model_name="model_GENE3D",
                   spm=spm,pswarm=pswarm,
                   mpi_ranks=1)
  with open("strike-slip.opts","w") as f:
    f.write(model.options)

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
  As an exemple for lithospheric models, the approximate Moho depth is a decent candidate.  

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
