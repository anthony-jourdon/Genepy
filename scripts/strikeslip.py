import numpy as np
import bcpy as bp

def material_parameters0():
  regions = {38 : [],
             39 : [],
             40 : [],
             41 : []}
  
  regions[38].append(bp.DensityBoussinesq("model_GENE3D",38,2700.0,3.0e-5,1.0e-11))
  regions[38].append(bp.PlasticDruckerPrager("model_GENE3D",38))
  regions[38].append(bp.SofteningLinear("model_GENE3D",38,0.0,0.5))
  regions[38].append(bp.ViscosityArrhenius2("model_GENE3D",38,"Quartzite"))

  regions[39].append(bp.DensityBoussinesq("model_GENE3D",39,2850.0,3.0e-5,1.0e-11))
  regions[39].append(bp.PlasticDruckerPrager("model_GENE3D",39))
  regions[39].append(bp.SofteningLinear("model_GENE3D",39,0.0,0.5))
  regions[39].append(bp.ViscosityArrhenius2("model_GENE3D",39,"Anorthite",Vmol=38.0e-6))

  regions[40].append(bp.DensityBoussinesq("model_GENE3D",40,3300.0,3.0e-5,1.0e-11))
  regions[40].append(bp.PlasticDruckerPrager("model_GENE3D",40))
  regions[40].append(bp.SofteningLinear("model_GENE3D",40,0.0,0.5))
  regions[40].append(bp.ViscosityArrhenius2("model_GENE3D",40,"Peridotite(dry)",Vmol=8.0e-6))

  regions[41].append(bp.DensityBoussinesq("model_GENE3D",41,3300.0,3.0e-5,1.0e-11))
  regions[41].append(bp.PlasticDruckerPrager("model_GENE3D",41))
  regions[41].append(bp.SofteningLinear("model_GENE3D",41,0.0,0.5))
  regions[41].append(bp.ViscosityArrhenius2("model_GENE3D",41,"Peridotite(dry)",Vmol=8.0e-6))

  opt = "########### Material parameters ###########\n"
  for region in regions:
    opt += f"###### Region {region} ######\n"
    for r in regions[region]:
      opt += r.sprint_option()
  return opt 

def material_parameters():
  regions = [
    bp.Region(38,
              bp.DensityBoussinesq(2700.0,3.0e-5,1.0e-11),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.ViscosityArrhenius2("Quartzite"),
              bp.Energy(1.5e-6,2.7)),
    bp.Region(39,
              bp.DensityBoussinesq(2850.0,3.0e-5,1.0e-11),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.ViscosityArrhenius2("Anorthite",Vmol=38.0e-6),
              bp.Energy(0.5e-6,2.85)),
    bp.Region(40,
              bp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              bp.Energy(0.0,3.3)),
    bp.Region(41,
              bp.DensityBoussinesq(3300.0,3.0e-5,1.0e-11),
              bp.SofteningLinear(0.0,0.5),
              bp.PlasticDruckerPrager(),
              bp.ViscosityArrhenius2("Peridotite(dry)",Vmol=8.0e-6),
              bp.Energy(0.0,3.3))
  ]

  all_regions = bp.ModelRegions("model_GENE3D",regions)
  opt = all_regions.sprint_option()
  return opt

def strikeslip():
  # 3D domain
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([16,16,16],      dtype=np.int32)   # Number of nodes
  # Create domain object
  Domain = bp.Domain(3,O,L,n)

  # Rotation of the referential
  r_angle = np.deg2rad(-15.0)                   # Rotation angle
  axis    = np.array([0,1,0], dtype=np.float64) # Rotation axis
  # Create rotation object
  Rotation = bp.Rotation(3,r_angle,axis)

  # velocity function
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0) # cm/a to m/s conversion
  u_norm  = 1.0 * cma2ms                   # horizontal velocity norm
  u_angle = np.deg2rad(90.0)               # velocity angle \in [-pi/2, pi/2]
  u_dir   = "z"                            # direction in which velocity varies
  u_type  = "extension"                    # extension or compression
  # Create boundary conditions object
  bc = bp.BoundaryConditions(Domain,u_norm,u_dir,u_type,u_angle,Rotation)

  # Evaluate the velocity and its derivatives
  u,grad_u = bc.evaluate_velocity_and_derivatives_symbolic() # symbolic
  u_num    = bc.evaluate_velocity_numeric()                  # numeric
  uL       = bc.get_velocity_orientation(horizontal=True,normalize=True)
  print(bc.report_symbolic_functions(u,grad_u,uL))

  # mesh refinement
  refinement = {"y": {"x_initial": np.array([-250,-180,-87.5,0], dtype=np.float64)*1e3,
                      "x_refined": np.array([-250,-50,-16.25,0], dtype=np.float64)*1e3}}
  m = bp.MeshRefinement(bc,refinement)
  m.refine()
  print(m)

  # gaussian initial strain
  ng = np.int32(2) # number of gaussians
  A  = np.array([1.0, 1.0],dtype=np.float64)
  # shape of the gaussians
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
  # Create gaussian object
  Gaussian = bp.Gaussian(m,Rotation,ng,A,a,b,c,x0,z0)
  Gaussian.evaluate_gaussians()
  print(Gaussian.report_symbolic_functions())
  strain = Gaussian.compute_field_distribution()

  # write the results to a file
  point_data = {"u": u_num, "strain": strain}
  w = bp.WriteVTS(m, vtk_fname="strike-slip.vts", point_data=point_data)
  w.write_vts()

  opt = Gaussian.sprint_option("model_GENE3D")
  opt += bc.sprint_option_dirichlet("model_GENE3D","Zfaces",43,["x","z"],u)
  opt += bc.sprint_option_navier("model_GENE3D","Xfaces",32,grad_u,uL)

  opt += material_parameters()
  print(opt)

if __name__ == "__main__":
  strikeslip()