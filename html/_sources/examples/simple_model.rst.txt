Example: simple model, linear viscous rheology
..............................................
This example shows how to build a very simple model using default linear viscous rheology 
and standard Dirichlet type boundary conditions for a velocity field imposing shortening in
:math:`z` direction as shown by the figure below.
Because the viscosity chosen in this example does not depend on temperature,
the thermal part of the model is not included.

.. image:: ../figures/compression.png
   :width: 200
   :align: center


1. Create a domain
~~~~~~~~~~~~~~~~~~~
We define a 3D :py:class:`domain <genepy.Domain>` :math:`\Omega = [0,600]\times[-250,0]\times[0,300]` km\ :sup:`3`
:math:`\in \mathbb R^3` discretized by a regular grid of 9x9x9 nodes.

.. code-block:: python

  import os
  import numpy as np
  import genepy as gp

  # 3D domain
  dimensions = 3
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = gp.Domain(dimensions,O,L,n)

2. Velocity function
~~~~~~~~~~~~~~~~~~~~
We define a simple orthogonal shortening :py:class:`velocity <genepy.VelocityLinear>` in the :math:`z` direction.

.. code-block:: python

  # velocity
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "compression"                  # extension or compression
  # Create Velocity class instance
  BCs = gp.VelocityLinear(Domain,u_norm,u_dir,u_type)

  # Access the symbolic velocity function
  u = BCs.u

.. note:: In this example, the derivatives of the velocity are not used.

3. Initial conditions
~~~~~~~~~~~~~~~~~~~~~
In this example we do not impose any initial plastic strain value nor mesh refinement.
Therefore the :py:class:`initial conditions <genepy.InitialConditions>` are only the Domain and the velocity function.
They will be used to generate the options for `pTatin3d`_ model.

.. code-block:: python

  # Initial conditions
  model_ics = gp.InitialConditions(Domain,u)

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
:doc:`boundary conditions <../boundary_conditions>` section.

.. code-block:: python

  # boundary conditions
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  u_bcs = [
    gp.Dirichlet(23,"Zmax",["z"],u, mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")), # orthogonal shortening
    gp.Dirichlet(37,"Zmin",["z"],u, mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")), # orthogonal shortening
    gp.Dirichlet(32,"Xmax",["x"],u, mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")), # free-slip
    gp.Dirichlet(14,"Xmin",["x"],u, mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")), # free-slip
    gp.DirichletUdotN(33,"Bottom",  mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")), # basal outflow
  ]
  # collect all boundary conditions
  model_bcs = gp.ModelBCs(u_bcs)

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
    gp.Region(38),
    # Lower crust
    gp.Region(39),
    # Lithosphere mantle
    gp.Region(40),
    # Asthenosphere
    gp.Region(41)
  ]
  model_regions = gp.ModelRegions(regions,
                                  mesh_file=os.path.join(root,"box_ptatin_md.bin"),
                                  region_file=os.path.join(root,"box_ptatin_region_cell.bin"))

6. Create the model and generate options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, we create the :py:class:`model <genepy.Model>` by gathering all the information defined previously and we save
the options to a file named ``simple_shortening_model.opts``.

.. code-block:: python

  # create class instance
  model = gp.Model(model_ics,model_regions,model_bcs)
  # write the options for ptatin3d
  with open("simple_shortening_model.opts","w") as f:
    f.write(model.options)
