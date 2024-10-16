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

def initial_strain_double_wz(Domain:gp.Domain,MshRef,Rotation,report=False):
  # gaussian initial strain
  ng = np.int32(2) # number of gaussians
  A  = np.array([1.0, 1.0],dtype=np.float64)
  # shape of the gaussians
  coeff = 0.5 * 6.0e-5**2
  a = np.array([coeff, coeff], dtype=np.float64)
  b = np.array([0.0, 0.0],     dtype=np.float64)
  c = np.array([coeff, coeff], dtype=np.float64)
  # position of the centre of the gaussians
  dz    = 27.5e3                         # distance from the domain centre in z direction
  angle = np.deg2rad(75)                 # angle between the x-axis and the line that passes through the centre of the domain and the centre of the gaussian
  domain_centre = 0.5*(Domain.O + Domain.L) # centre of the domain
  
  x0 = np.zeros(shape=(ng), dtype=np.float64)
  # centre of the gaussian in z direction
  z0 = np.array([domain_centre[2] - dz, 
                 domain_centre[2] + dz], dtype=np.float64) 
  # centre of the gaussian in x direction
  x0[0] = gp.utils.x_centre_from_angle(z0[0],angle,(domain_centre[0],domain_centre[2])) 
  x0[1] = gp.utils.x_centre_from_angle(z0[1],angle,(domain_centre[0],domain_centre[2]))
  # Create gaussian object
  Gaussian = gp.Gaussian(MshRef,ng,A,a,b,c,x0,z0,Rotation)
  Gaussian.evaluate_gaussians()
  if report:
    print(Gaussian.report_symbolic_functions())
  strain = Gaussian.compute_field_distribution()
  return Gaussian,strain

def initial_strain_single_wz(Domain:gp.Domain,MshRef,Rotation,report=False):
  # gaussian initial strain
  ng = np.int32(1) # number of gaussians
  A  = np.array([1.0],dtype=np.float64)
  # shape of the gaussians
  coeff = 0.5 * 6.0e-5**2
  a = np.array([coeff], dtype=np.float64)
  b = np.array([0.0],     dtype=np.float64)
  c = np.array([coeff], dtype=np.float64)
  # position of the centre of the gaussian
  domain_centre = 0.5*(Domain.O + Domain.L) # centre of the domain
  x0 = np.array([domain_centre[0]], dtype=np.float64) # centre of the gaussian in x direction
  z0 = np.array([domain_centre[2]], dtype=np.float64) # centre of the gaussian in z direction
  # Create gaussian object
  Gaussian = gp.Gaussian(MshRef,ng,A,a,b,c,x0,z0,Rotation)
  Gaussian.evaluate_gaussians()
  if report:
    print(Gaussian.report_symbolic_functions())
  strain = Gaussian.compute_field_distribution()
  return Gaussian,strain

def initial_strain_multiple_gaussians_aligned(Domain:gp.Domain,MshRef,Rotation,report=False):
  # gaussian initial strain
  ng = np.int32(2) # number of gaussians
  A  = np.array([1.0, 1.0],dtype=np.float64)
  # shape of the gaussians
  coeff = 0.5 * 6.0e-5**2
  a = np.array([coeff, coeff], dtype=np.float64)
  b = np.array([0.0, 0.0],     dtype=np.float64)
  c = np.array([coeff, coeff], dtype=np.float64)

  domain_centre = 0.5*(Domain.O_num + Domain.L_num) # centre of the domain

  z0 = np.array([domain_centre[2], domain_centre[2]], dtype=np.float64)
  x0 = np.array([100.0e3, 500.0e3], dtype=np.float64)
  
  Gaussian = gp.Gaussian(MshRef,ng,A,a,b,c,x0,z0,Rotation)
  strain = Gaussian.compute_field_distribution()
  return Gaussian,strain

def initial_heat_source(Domain,Rotation,report=False):
  # gaussian initial heat source
  dc = 50.0e3 # distance between gaussian centres
  # number of gaussians
  nx = int(Domain.L[0] / dc)
  nz = int(Domain.L[2] / dc)

  xc = np.zeros(shape=(nx*nz),dtype=np.float64)
  zc = np.zeros(shape=(nx*nz),dtype=np.float64)
  n = 0
  for j in range(nz):
    for i in range(nx):
      xc[n] = dc * (i + np.random.rand())
      zc[n] = dc * (j + np.random.rand())
      n += 1

  # shape of the gaussians
  coeff = 0.5 * 6.0e-5**2
  a = np.ones(shape=(nx*nz),dtype=np.float64) * coeff
  b = np.zeros(shape=(nx*nz),dtype=np.float64)
  c = np.ones(shape=(nx*nz),dtype=np.float64) * coeff

  # amplitude of the gaussians (heat source)
  A = np.ones(shape=(nx*nz),dtype=np.float64) * 1.5e-6

  # Create gaussian object
  Gaussian = gp.Gaussian(Domain,nx*nz,A,a,b,c,xc,zc)
  Gaussian.evaluate_gaussians()
  return Gaussian

def initial_conditions(Domain,MshRef,IniStrain,u):
  plstr = gp.InitialPlasticStrain(IniStrain)
  model_ics = gp.InitialConditions(Domain,u,mesh_refinement=MshRef,initial_strain=plstr)
  return model_ics

def boundary_conditions(u,grad_u,uL):
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  #root = os.path.join(os.environ['SOFTS'],"gmsh_to_point_cloud")
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

def boundary_conditions_navier_slip_all(u,grad_u,uL):
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  bcs = [
    gp.Dirichlet(tag=23,name="Zmax",components=["x","z"],velocity=u,mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")),
    gp.NavierSlip(37,"Zmin",grad_u,uL,mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")),
    gp.NavierSlip(tag=32,name="Xmax",grad_u=grad_u,u_orientation=uL,mesh_file=os.path.join(root,"box_ptatin_facet_32_mesh.bin")),
    gp.NavierSlip(14,"Xmin",grad_u,uL,mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
    gp.DirichletUdotN(33,"Bottom",mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")),
  ]
  # Temperature boundary conditions
  Tbcs = gp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  all_bcs = gp.ModelBCs(bcs,Tbcs)
  return all_bcs

def boundary_conditions_compose(u,grad_u,uL):
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  bcs = [
    gp.Dirichlet(tag=23,name="Zmax",components=["x","z"],velocity=u,mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")),
    gp.Dirichlet(37,"Zmin",["x","z"],u,mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")),
    gp.NavierSlip(tag=32,name="Xmax",grad_u=grad_u,u_orientation=uL,mesh_file=os.path.join(root,"box_ptatin_facet_32_mesh.bin")),
    gp.NavierSlip(14,"Xmin",grad_u,uL,mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
    gp.Dirichlet(tag=12,name="XminDiri",components=["y"],velocity=u,mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
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
              gp.Energy(heat_source=gp.EnergySource(#gp.EnergySourceMaterialPointValue(),
                                                    gp.EnergySourceConstant(1.5e-6),
                                                    gp.EnergySourceShearHeating()),
                        conductivity=2.7)),
    # Lower crust
    gp.Region(39,
              gp.DensityBoussinesq(density=2850.0,thermal_expansion=3.0e-5,compressibility=1.0e-11),
              gp.ViscosityArrhenius2("Quartzite"),
              gp.SofteningLinear(strain_min=0.0,strain_max=0.5),
              gp.PlasticDruckerPrager(),
              gp.Energy(heat_source=gp.EnergySource(gp.EnergySourceConstant(0.5e-6),
                                                    gp.EnergySourceShearHeating()),
                        conductivity=2.85)),
    # Lithosphere mantle
    gp.Region(40,
              gp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              gp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              gp.SofteningLinear(0.0,0.5),
              gp.PlasticDruckerPrager(),
              gp.Energy(heat_source=gp.EnergySource(gp.EnergySourceShearHeating()),
                        conductivity=3.3)),
    # Asthenosphere
    gp.Region(41,
              gp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              gp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              gp.SofteningLinear(0.0,0.5),
              gp.PlasticDruckerPrager(),
              gp.Energy(heat_source=gp.EnergySource(gp.EnergySourceShearHeating()),
                        conductivity=3.3))
  ]

  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  #root = os.path.join(os.environ['SOFTS'],"gmsh_to_point_cloud")
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
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  all_regions = gp.ModelRegions(regions,
                                mesh_file=os.path.join(root,"box_ptatin_md.bin"),
                                region_file=os.path.join(root,"box_ptatin_region_cell.bin"))
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
  Gaussian,strain = initial_strain_double_wz(Domain,MshRef,Rotation,report=True)
  #Gaussian,strain = initial_strain_single_wz(Domain,MshRef,Rotation,report=True)
  #Gaussian,strain = initial_strain_multiple_gaussians_aligned(Domain,MshRef,Rotation,report=True)
  # initial heat source
  Gaussian_hs = initial_heat_source(Domain,Rotation,report=False)

  #plstr  = gp.InitialPlasticStrain(Gaussian)
  #hs_ini = gp.InitialHeatSource(Gaussian_hs)
  #ics    = gp.InitialConditions(Domain,BCs.u,mesh_refinement=MshRef,initial_strain=plstr,initial_heat_source=hs_ini)

  # write the results to a vts file for visualization
  #point_data = {"u": u_num, "strain": strain}
  #w = gp.WriteVTS(MshRef, vtk_fname="strike-slip.vts", point_data=point_data)
  #w.write_vts()

  # generate objects for options writing
  ics     = initial_conditions(Domain,MshRef,Gaussian,BCs.u)
  bcs     = boundary_conditions(BCs.u,BCs.grad_u,BCs.u_dir_horizontal)
  #bcs     = boundary_conditions_compose(BCs.u,BCs.grad_u,BCs.u_dir_horizontal)
  #bcs     = boundary_conditions_navier_slip_all(BCs.u,BCs.grad_u,BCs.u_dir_horizontal)
  regions = material_parameters()
  #regions = test_default_material_parameters()
  spm = gp.SPMDiffusion(["zmin","zmax"],diffusivity=1.0e-6)
  pswarm = gp.PswarmFillBox([0.0,-100.0e3,0.0],[600e3,-4.0e3,300.0e3],layout=[30,5,15],pressure=True,temperature=True)
  markers = gp.MarkersManagement(layout=(8,8,8),
                                 popctrl_faces=(0, 1, 4, 5), 
                                 popctrl_np_lower=8,
                                 popctrl_np_upper=128, 
                                 popctrl_layout=(2, 2, 2))

  # write the options for ptatin3d
  model = gp.Model(ics,regions,bcs,
                   model_name="model_GENE3D",
                   spm=spm,#pswarm=pswarm,
                   markers=markers,
                   output_fields=["region","viscosity","density","plastic_strain"],
                   mpi_ranks=1)
  #print(model.options)
  with open("strike-slip.sh","w") as f:
    f.write(model.options)

if __name__ == "__main__":
  strikeslip()