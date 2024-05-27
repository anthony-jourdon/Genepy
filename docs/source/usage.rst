.. _pTatin3d: https://github.com/laetitialp/ptatin-gene

Usage
=====

.. _installation:

Installation
------------
First, this module requires the following dependencies that can all be installed using pip

- numpy

- sympy

- matplotlib

- pyvista

To be able to import bcpy module its location needs to be known by 
the environment variable **$PYTHONPATH**. 
To do so, move yourself in bc-pre-processing and type

.. code-block:: console

  source bc-pre-proc.sh

This needs to be done every time a new terminal or python session is opened.
Another possibility is to add it to your bashrc or zshrc but it may be erased by python virtual environment when loaded.

Once your **$PYTHONPATH** has been appended check that it works correctly by typing

.. code-block:: python

  import bcpy

in a python terminal. If no error is raised, the module is correctly installed.

Building a model
----------------

Exemple: strike-slip model
..........................
This example will build a 3D model with vertical mesh refinement and a strike-slip velocity
field rotated by 15 degrees as showed in the figure below.

.. image:: figures/Strike_slip-01.png
   :width: 600
   :align: center

.. note:: 
  The following code blocks may re-use variables defined in previous blocks.
  To ensure that the code runs correctly, 
  it is recommended to run the code blocks in the order they are presented.

1. Create a domain
~~~~~~~~~~~~~~~~~~~
We define a 3D domain :math:`\Omega = [0,600]\times[-250,0]\times[0,300]` km\ :sup:`3` 
:math:`\in \mathbb R^3` discretized by a regular grid of 9x9x9 nodes. 

.. code-block:: python

  import numpy as np
  import bcpy as bp

  # 3D domain
  dimensions = 3
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = bp.Domain(dimensions,O,L,n)

2. Mesh refinement
~~~~~~~~~~~~~~~~~~
This step is optional. 
It allows to refine the mesh in any direction using linear interpolation.
In this example we refine the mesh in the y direction.

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
To rotate the velocity field we first need to set the parameters of this rotation.
This is achieved with the Rotation class.

.. code-block:: python

  # Rotation of the referential
  r_angle = np.deg2rad(-15.0)                   # Rotation angle
  axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
  # Create instance of Rotation class
  Rotation = bp.Rotation(dimensions,r_angle,axis)

4. Velocity field
~~~~~~~~~~~~~~~~~
The following code block creates a strike-slip velocity field with a norm of 1 cm/a.
The method 
:py:meth:`evaluate_velocity_and_derivatives_symbolic() <bcpy.Velocity.evaluate_velocity_and_derivatives_symbolic>` 
returns the symbolic expression of the velocity field and its gradient.
The method
:py:meth:`evaluate_velocity_numeric() <bcpy.Velocity.evaluate_velocity_numeric>`
returns the numeric value of the velocity field evaluated at coordinates of the nodes.
The method
:py:meth:`get_velocity_orientation() <bcpy.Velocity.get_velocity_orientation>`
returns the orientation of the velocity field at the boundary.

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
  u,grad_u = BCs.evaluate_velocity_and_derivatives_symbolic() # symbolic
  u_num    = BCs.evaluate_velocity_numeric()                  # numeric
  # Get the orientation of the vectors at boundary (horizontal removes the vertical component)
  uL       = BCs.get_velocity_orientation(horizontal=True,normalize=True)

5. Define gaussian weak zones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this exemple we define two gaussian weak zones.
This step is not mandatory if your model already contains a 
way to localize deformation.

.. code-block:: python

  # gaussian weak zones
  ng = np.int32(2) # number of gaussians
  A  = np.array([1.0, 1.0],dtype=np.float64) # amplitude (will be multiplied by a random number to generate noise in the model)
  # coefficients of the shape of the gaussians
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
In the following example a path to the mesh files describing the boundaries is provided.
These mesh files are located in ``"ptatin-gene/src/models/gene3d/examples"``.
You can modify the ``root`` variable to match the location of the mesh files on your system 
or remove that part of the code if you do not have access to these files.
Note however that `pTatin3d`_ requires mesh files to define the boundaries.

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
  Tbcs = bp.TemperatureBC(faces=["ymax","ymin"],values=[0.0,1450.0])
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

9. Create the model and generate options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The model is created by gathering all the information defined previously plus some 
additional parameters.

.. code-block:: python

  # Add erosion-sedimentation with diffusion
  spm = bp.SPMDiffusion(["zmin","zmax"],diffusivity=1.0e-6)
  # Add passive tracers
  pswarm = bp.PswarmFillBox([0.0,-100.0e3,0.0],[600e3,-4.0e3,300.0e3],layout=[30,5,15],pressure=True,temperature=True)

  # write the options for ptatin3d
  model = bp.Model(model_ics,model_regions,model_bcs,
                   model_name="model_GENE3D",
                   spm=spm,pswarm=pswarm,
                   mpi_ranks=1)
  with open("strike-slip.opts","w") as f:
    f.write(model.options)

