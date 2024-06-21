import os
import numpy as np
import genepy as gp

def model_domain():
  # 3D domain
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create domain object
  Domain = gp.Domain(3,O,L,n)
  return Domain

def mesh_refinement(Domain,report=False):
  # mesh refinement
  refinement = {"y": {"x_initial": np.array([-250,-180,-87.5,0], dtype=np.float64)*1e3,
                      "x_refined": np.array([-250,-50,-16.25,0], dtype=np.float64)*1e3}}
  MshRef = gp.MeshRefinement(Domain,refinement)
  MshRef.refine()
  if report:
    print(MshRef)
  return MshRef

def domain_rotation():
  # Rotation of the referential
  r_angle = np.deg2rad(-15.0)                   # Rotation angle
  axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
  # Create rotation object
  Rotation = gp.Rotation(3,r_angle,axis)
  return Rotation

def velocity_bcs(Domain,Rotation,report=False):
  # velocity function
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_angle = np.deg2rad(90.0)                # velocity angle \in [-pi/2, pi/2]
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "extension"                    # extension or compression
  # Create boundary conditions object
  BCs = gp.VelocityLinear(Domain,u_norm,u_dir,u_type,u_angle,Rotation)

  # Evaluate the velocity and its derivatives
  u_num    = BCs.evaluate_velocity_numeric() # numeric
  if report:
    print(BCs.report_symbolic_functions())
  return BCs,u_num

def initial_strain(Domain,MshRef,Rotation,report=False):
  # gaussian initial strain
  ng = np.int32(2) # number of gaussians
  A  = np.array([1.0, 1.0],dtype=np.float64)
  # shape of the gaussians
  coeff = 0.5 * 6.0e-5**2
  a = np.array([coeff, coeff], dtype=np.float64)
  b = np.array([0.0, 0.0],     dtype=np.float64)
  c = np.array([coeff, coeff], dtype=np.float64)
  # position of the centre of the gaussians
  dz    = 30.0e3                            # distance from the domain centre in z direction
  angle = np.deg2rad(-70.0)                 # angle between the x-axis and the line that passes through the centre of the domain and the centre of the gaussian
  domain_centre = 0.5*(Domain.O + Domain.L) # centre of the domain
  
  x0 = np.zeros(shape=(ng), dtype=np.float64)
  # centre of the gaussian in z direction
  z0 = np.array([domain_centre[2] - dz, 
                 domain_centre[2] + dz], dtype=np.float64) 
  # centre of the gaussian in x direction
  x0[0] = gp.utils.x_centre_from_angle(z0[0],angle,(domain_centre[0],domain_centre[2])) 
  x0[1] = gp.utils.x_centre_from_angle(z0[1],angle,(domain_centre[0],domain_centre[2]))
  # Create gaussian object
  Gaussian = gp.Gaussian(MshRef,Rotation,ng,A,a,b,c,x0,z0)
  Gaussian.evaluate_gaussians()
  if report:
    print(Gaussian.report_symbolic_functions())
  strain = Gaussian.compute_field_distribution()
  return Gaussian,strain

def initial_conditions(Domain,MshRef,IniStrain,u):
  model_ics = gp.InitialConditions(Domain,u,mesh_refinement=MshRef,initial_strain=IniStrain)
  return model_ics

def boundary_conditions(u,grad_u,uL):
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  bcs = [
    gp.Dirichlet(tag=23,name="Zmax",components=["x","z"],velocity=u,mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")),
    gp.Dirichlet(37,"Zmin",["x","z"],u,mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")),
    gp.NavierSlip(tag=32,name="Xmax",grad_u=grad_u,u_orientation=uL,mesh_file=os.path.join(root,"box_ptatin_facet_32_mesh.bin")),
    gp.NavierSlip(14,"Xmin",grad_u,uL,mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
    gp.DirichletUdotN(33,"Bottom",mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")),
  ]
  # Temperature boundary conditions
  Tbcs = gp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  all_bcs = gp.ModelBCs(bcs,Tbcs)
  return all_bcs

def material_parameters():
  # Define the material parameters for the model as a list of Region objects
  regions = [
    # Upper crust
    gp.Region(38,                                          # region tag
              gp.DensityBoussinesq(2700.0,3.0e-5,1.0e-11), # density
              gp.ViscosityArrhenius2("Quartzite"),         # viscosity  (values from the database using rock name)
              gp.SofteningLinear(0.0,0.5),                 # softening
              gp.PlasticDruckerPrager(),                   # plasticity (default values, can be modified using the corresponding parameters)
              gp.Energy(1.5e-6,2.7)),                      # energy
    # Lower crust
    gp.Region(39,
              gp.DensityBoussinesq(density=2850.0,thermal_expansion=3.0e-5,compressibility=1.0e-11),
              gp.ViscosityArrhenius2("Anorthite",Vmol=38.0e-6),
              gp.SofteningLinear(strain_min=0.0,strain_max=0.5),
              gp.PlasticDruckerPrager(),
              gp.Energy(heat_source=0.5e-6,conductivity=2.85)),
    # Lithosphere mantle
    gp.Region(40,
              gp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              gp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              gp.SofteningLinear(0.0,0.5),
              gp.PlasticDruckerPrager(),
              gp.Energy(0.0,3.3)),
    # Asthenosphere
    gp.Region(41,
              gp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              gp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              gp.SofteningLinear(0.0,0.5),
              gp.PlasticDruckerPrager(),
              gp.Energy(0.0,3.3))
  ]

  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  all_regions = gp.ModelRegions(regions,
                                mesh_file=os.path.join(root,"box_ptatin_md.bin"),
                                region_file=os.path.join(root,"box_ptatin_region_cell.bin"))
  return all_regions

def test_default_material_parameters():
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
  all_regions = gp.ModelRegions(regions)
  return all_regions

def strikeslip():
  # model domain
  Domain   = model_domain()
  # rotation of the referential, if not needed, the Velocity object can be created without it
  Rotation = domain_rotation()
  # boundary conditions
  BCs,u_num = velocity_bcs(Domain,Rotation,report=True)
  # mesh refinement
  MshRef = mesh_refinement(BCs,report=False)
  # initial strain
  Gaussian,strain = initial_strain(Domain,MshRef,Rotation,report=True)

  # write the results to a vts file for visualization
  point_data = {"u": u_num, "strain": strain}
  w = gp.WriteVTS(MshRef, vtk_fname="strike-slip.vts", point_data=point_data)
  w.write_vts()

  # generate objects for options writing
  ics     = initial_conditions(Domain,MshRef,Gaussian,BCs.u)
  bcs     = boundary_conditions(BCs.u,BCs.grad_u,BCs.u_dir_horizontal)
  regions = material_parameters()
  #regions = test_default_material_parameters()
  spm = gp.SPMDiffusion(["zmin","zmax"],diffusivity=1.0e-6)
  pswarm = gp.PswarmFillBox([0.0,-100.0e3,0.0],[600e3,-4.0e3,300.0e3],layout=[30,5,15],pressure=True,temperature=True)

  # write the options for ptatin3d
  model = gp.Model(ics,regions,bcs,
                   model_name="model_GENE3D",
                   spm=spm,pswarm=pswarm,
                   mpi_ranks=1)
  #print(model.options)
  with open("strike-slip.sh","w") as f:
    f.write(model.options)

if __name__ == "__main__":
  strikeslip()