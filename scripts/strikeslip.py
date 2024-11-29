import os
import numpy as np
import genepy as gp
import matplotlib.pyplot as plt

def plot_gaussian(Domain:gp.Domain,gaussians:list):
  field = np.zeros(shape=(Domain.num_coor[0].shape), dtype=np.float64)
  for g in gaussians:
    g.evaluate_gaussian()
    field += g.gaussian_num

  _,ax = plt.subplots(tight_layout=True)
  c = ax.contourf(Domain.num_coor[0][:,0,:],Domain.num_coor[2][:,0,:],field[:,0,:], levels=100, cmap='plasma')
  ax.set_xlim([Domain.L_num[0],Domain.O_num[0]])
  ax.axis('equal')
  ax.set_title('Field distribution')
  ax.set_xlabel('x')
  ax.set_ylabel('z')
  plt.colorbar(c, ax=ax, orientation='horizontal')
  plt.draw()
  return

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

class InitialWeakZones:
  def __init__(
      self, 
      type:str, 
      Domain:gp.Domain, 
      coeff:float,
      dz:float=0.0,
      angle:float=0.0,
      Rotation:gp.Rotation|None=None,
    ):
    self.type     = type
    self.Domain   = Domain
    self.Rotation = Rotation
    self.coeff    = coeff
    self.dz       = dz
    self.angle    = angle

    self.centre = 0.5*(self.Domain.O_num + self.Domain.L_num) # centre of the domain
    self.Gaussian = None
    self.create_gaussian()
    return
  
  def create_gaussian(self):
    if self.type == "single":
      self.single_gaussian()
    elif self.type == "double":
      self.double_gaussian()
    elif self.type == "triple":
      self.triple_gaussians()
    else:
      raise ValueError("Invalid type of weak zone, choose between 'single', 'double' or 'triple', found: {}".format(self.type))
    return

  def single_gaussian(self):
    # Create gaussian object
    Gaussian = gp.Gaussian2D(self.Domain,1.0,self.coeff,self.coeff,self.centre[0],self.centre[2],self.Rotation)
    self.Gaussian = [Gaussian]
    return 

  def double_gaussian(self):
    ng = np.int32(2) # number of gaussians
    
    x0 = np.zeros(shape=(ng), dtype=np.float64)
    # centre of the gaussian in z direction
    z0 = np.array(
      [
        self.centre[2] - self.dz, 
        self.centre[2] + self.dz
      ], dtype=np.float64
    ) 
    # centre of the gaussian in x direction
    for i in range(ng):
      x0[i] = gp.utils.x_centre_from_angle(z0[i],self.angle,(self.centre[0],self.centre[2]))
    # Create gaussian object
    Gaussian = []
    for i in range(ng):
      Gaussian.append(gp.Gaussian2D(self.Domain,1.0,self.coeff,self.coeff,x0[i],z0[i],self.Rotation))
    self.Gaussian = Gaussian
    return

  def triple_gaussians(self):
    # gaussian initial strain
    ng = np.int32(3) # number of gaussians
    x0 = np.zeros(shape=(ng), dtype=np.float64)
    # centre of the gaussian in z direction
    z0 = np.array(
      [
        self.centre[2] - self.dz, 
        self.centre[2] + self.dz,
        self.centre[2]
      ], dtype=np.float64
    ) 
    # centre of the gaussian in x direction
    for i in range(ng-1):
      x0[i] = gp.utils.x_centre_from_angle(z0[i],self.angle,(self.centre[0],self.centre[2]))
    x0[2] = self.centre[0]
    # Create gaussian object
    Gaussian = []
    for i in range(ng):
      Gaussian.append(gp.Gaussian2D(self.Domain,1.0,self.coeff,self.coeff,x0[i],z0[i],self.Rotation))
    self.Gaussian = Gaussian
    return

def model_domain(
    Ox:float, Oy:float, Oz:float,
    Lx:float, Ly:float, Lz:float,
    nx:int,   ny:int,   nz:int
  ) -> gp.Domain:
  # 3D domain
  O = np.array([Ox, Oy, Oz], dtype=np.float64) # Origin
  L = np.array([Lx, Ly, Lz], dtype=np.float64) # Length
  n = np.array([nx, ny, nz], dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create domain object
  Domain = gp.Domain(3,O,L,n)
  return Domain

def main(report=False,plot=False):	
  # model Domain
  Domain = model_domain(
    Ox=0.0,   Oy=-250e3, Oz=0.0, 
    Lx=600e3, Ly=0.0,    Lz=300e3, 
    nx=128,   ny=64,     nz=256
  )
  
  # Rotation of the referential
  r_angle = np.deg2rad(-15.0)                   # Rotation angle
  axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
  # Create rotation object
  Rotation = gp.Rotation(3,r_angle,axis)
  
  # velocity function
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_angle = np.deg2rad(90.0)                # velocity angle \in [-pi/2, pi/2]
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "extension"                    # extension or compression
  # Create boundary conditions object
  BCs = gp.VelocityLinear(Domain,u_norm,u_dir,u_type,u_angle,Rotation)
  if report:
    print(BCs.report_symbolic_functions())
  
  # mesh refinement
  refinement = {"y": {"x_initial": np.array([-250,-180,-87.5,0], dtype=np.float64)*1e3,
                      "x_refined": np.array([-250,-50,-16.25,0], dtype=np.float64)*1e3}}
  Mesh_refinement = gp.MeshRefinement(Domain,refinement)
  if report:
    print(Mesh_refinement)

  # initial strain
  coeff = 0.5*6.0e-5**2 # shape of the gaussians
  dz    = 25e3 # distance from the domain centre in z direction
  angle = np.deg2rad(0) # angle between the x-axis and the line that passes through the centre of the domain and the centre of the gaussian
  IS = InitialWeakZones("triple",Domain,coeff,dz,angle,Rotation)
  if report:
    for g in IS.Gaussian:
      print(g)
  if plot:
    plot_gaussian(Domain,IS.Gaussian)
  Gopts = gp.GaussiansOptions(IS.Gaussian)

  plstr = gp.InitialPlasticStrain(Gopts)
  model_ics = gp.InitialConditions(Domain,BCs.u,mesh_refinement=Mesh_refinement,initial_strain=plstr)

  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  bcs = [
    gp.Dirichlet(  tag=23, name="Zmax", components=["x","z"], velocity=BCs.u, mesh_file=os.path.join(root,"box_ptatin_facet_23_mesh.bin")),
    gp.Dirichlet(  tag=37, name="Zmin", components=["x","z"], velocity=BCs.u, mesh_file=os.path.join(root,"box_ptatin_facet_37_mesh.bin")),
    gp.NavierSlip( tag=32, name="Xmax", grad_u=BCs.grad_u, u_orientation=BCs.u_dir_horizontal, mesh_file=os.path.join(root,"box_ptatin_facet_32_mesh.bin")),
    gp.NavierSlip( tag=14, name="Xmin", grad_u=BCs.grad_u, u_orientation=BCs.u_dir_horizontal, mesh_file=os.path.join(root,"box_ptatin_facet_14_mesh.bin")),
    gp.DirichletUdotN(33,"Bottom",mesh_file=os.path.join(root,"box_ptatin_facet_33_mesh.bin")),
  ]
  # Temperature boundary conditions
  Tbcs = gp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  model_bcs = gp.ModelBCs(bcs,Tbcs)

  regions = material_parameters()
  spm = gp.SPMDiffusion(["zmin","zmax"],diffusivity=1.0e-6)
  pswarm = gp.PswarmFillBox([0.0,-100.0e3,0.0],[600e3,-4.0e3,300.0e3],layout=[30,5,15],pressure=True,temperature=True)
  markers = gp.MarkersManagement(layout=(16,16,16),
                                 popctrl_faces=(0, 1, 4, 5), 
                                 popctrl_np_lower=8,
                                 popctrl_np_upper=128, 
                                 popctrl_layout=(2, 2, 2))

  # write the options for ptatin3d
  model = gp.Model(model_ics,regions,model_bcs,
                   model_name="model_GENE3D",
                   spm=spm,#pswarm=pswarm,
                   markers=markers,
                   output_fields=["region","viscosity","density","plastic_strain","damage"],
                   mpi_ranks=1)
  #print(model.options)
  with open("strike-slip.sh","w") as f:
    f.write(model.options)
  return

if __name__ == "__main__":
  main(report=True,plot=False)
  plt.show()