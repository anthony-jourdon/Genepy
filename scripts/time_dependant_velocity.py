import os
import numpy as np
import genepy as gp

def extension(Domain,cma2ms):
  u_params  = {"u_norm"        : 1.0 * cma2ms,     # norm of the velocity on boundary
               "variation_dir" : "x",              # direction in which velocity varies 
               "velocity_type" : "extension",      # extension or compression (because norm > 0)
               "u_angle"       : np.deg2rad(90.0)  # angle of the velocity vector
              }
  # Create linear velocity class instance
  u_ext     = gp.VelocityLinear(Domain,**u_params)
  return u_ext

def compression(Domain,cma2ms):
  u_params  = {"u_norm"        : 3.0 * cma2ms,     # norm of the velocity on boundary
               "variation_dir" : "x",              # direction in which velocity varies 
               "velocity_type" : "compression",    # extension or compression (because norm > 0)
               "u_angle"       : np.deg2rad(45.0)  # angle of the velocity vector
              }
  # Create linear velocity class instance
  u_comp    = gp.VelocityLinear(Domain,**u_params)
  return u_comp

def boundary_conditions_Navier(u,grad_u,uL):
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  # zmax parameters
  zmax = {"tag"        :23, 
          "name"       :"Zmax", 
          "components" :["x","z"], 
          "velocity"   :u, 
          "mesh_file"  :os.path.join(root,"box_ptatin_facet_23_mesh.bin")}
  # Create Dirichlet boundary condition class instance for zmax face
  zmax_face = gp.Dirichlet(**zmax)

  # zmin parameters
  zmin = {"tag":37, 
          "name":"Zmin", 
          "components":["x","z"], 
          "velocity":u, 
          "mesh_file":os.path.join(root,"box_ptatin_facet_37_mesh.bin")}
  # Create Dirichlet boundary condition class instance for zmin face
  zmin_face = gp.Dirichlet(**zmin)
  
  # xmax parameters
  xmax = {"tag":32, 
          "name":"Xmax", 
          "grad_u":grad_u, 
          "u_orientation":uL, 
          "mesh_file":os.path.join(root,"box_ptatin_facet_32_mesh.bin")}
  # Create Navier slip boundary condition class instance for xmax face
  xmax_face = gp.NavierSlip(**xmax)
  
  # xmin parameters
  xmin = {"tag":14, 
          "name":"Xmin", 
          "grad_u":grad_u, 
          "u_orientation":uL, 
          "mesh_file":os.path.join(root,"box_ptatin_facet_14_mesh.bin")}
  # Create Navier slip boundary condition class instance for xmin face
  xmin_face = gp.NavierSlip(**xmin)
  
  # bottom parameters
  bottom = {"tag":33, 
            "name":"Bottom", 
            "mesh_file":os.path.join(root,"box_ptatin_facet_33_mesh.bin")}
  # Create DirichletUdotN boundary condition class instance for bottom face
  bottom_face = gp.DirichletUdotN(**bottom)

  # collect stokes boundary conditions in a list
  bcs = [zmax_face,zmin_face,xmax_face,xmin_face,bottom_face]
  # Temperature boundary conditions
  Tbcs = gp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  all_bcs = gp.ModelBCs(bcs,Tbcs)
  return all_bcs

def boundary_conditions_Dirichlet(u):
  # path to mesh files (system dependent, change accordingly)
  root = os.path.join(os.environ['PTATIN'],"ptatin-gene/src/models/gene3d/examples")
  # Velocity boundary conditions
  # zmax parameters
  zmax = {"tag"        :23, 
          "name"       :"Zmax", 
          "components" :["z"], 
          "velocity"   :u, 
          "mesh_file"  :os.path.join(root,"box_ptatin_facet_23_mesh.bin")}
  # Create Dirichlet boundary condition class instance for zmax face
  zmax_face = gp.Dirichlet(**zmax)

  # zmin parameters
  zmin = {"tag":37, 
          "name":"Zmin", 
          "components":["z"], 
          "velocity":u, 
          "mesh_file":os.path.join(root,"box_ptatin_facet_37_mesh.bin")}
  # Create Dirichlet boundary condition class instance for zmin face
  zmin_face = gp.Dirichlet(**zmin)
  
  # xmax parameters
  xmax = {"tag":32, 
          "name":"Xmax", 
          "components":["x"], 
          "velocity":u, 
          "mesh_file":os.path.join(root,"box_ptatin_facet_32_mesh.bin")}
  # Create Navier slip boundary condition class instance for xmax face
  xmax_face = gp.Dirichlet(**xmax)
  
  # xmin parameters
  xmin = {"tag":14, 
          "name":"Xmin", 
          "components":["x"], 
          "velocity":u, 
          "mesh_file":os.path.join(root,"box_ptatin_facet_14_mesh.bin")}
  # Create Navier slip boundary condition class instance for xmin face
  xmin_face = gp.Dirichlet(**xmin)
  
  # bottom parameters
  bottom = {"tag":33, 
            "name":"Bottom", 
            "mesh_file":os.path.join(root,"box_ptatin_facet_33_mesh.bin")}
  # Create DirichletUdotN boundary condition class instance for bottom face
  bottom_face = gp.DirichletUdotN(**bottom)

  # collect stokes boundary conditions in a list
  bcs = [zmax_face,zmin_face,xmax_face,xmin_face,bottom_face]
  # Temperature boundary conditions
  Tbcs = gp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
  # collect all boundary conditions
  all_bcs = gp.ModelBCs(bcs,Tbcs)
  return all_bcs

def main():
  # 3D domain
  dimensions = 3
  O = np.array([0,-250e3,0],    dtype=np.float64) # Origin
  L = np.array([600e3,0,300e3], dtype=np.float64) # Length
  n = np.array([9,9,9],         dtype=np.int32)   # Number of Q1 nodes i.e. elements + 1
  # Create Domain class instance
  Domain = gp.Domain(dimensions,O,L,n)

  # cm/a to m/s conversion
  cma2ms  = 1e-2 / (3600.0 * 24.0 * 365.0)

  phase_1 = extension(Domain,cma2ms)
  phase_2 = compression(Domain,cma2ms)
  
  # time inversion parameters
  Ma2s = (3600.0 * 24.0 * 365.0) * 1e6
  t1   = 2.0 * Ma2s
  t2   = 4.0 * Ma2s
  breakpoints = np.array([ t1, t2 ], dtype=np.float64)  # breakpoints in time (where atan(t-t0)=0 )
  slopes = np.array([ 5e-13, 5e-13 ], dtype=np.float64) # rate of change of atan(s(t-t0)) function

  # create time dependant velocity inversion class instance
  bc_inv = gp.VelocityInversion(Domain,phase_1,phase_2,breakpoints,slopes)
  # space and time dependant velocity function and its gradient
  u,grad_u = bc_inv.evaluate_velocity_and_gradient_symbolic()
  print(u,grad_u)
  # Orientation of the velocity at the boundary 
  uL_1,uL_2 = bc_inv.get_velocity_orientation(horizontal=True,normalize=True)
  # Report the computed function and its gradient
  print(bc_inv.report_symbolic_functions(u,grad_u))

  # time at which velocity is 0 during the tectonic regime inversion
  t0 = bc_inv.get_time_zero_velocity(report=True)
  print(f"Time at which velocity is 0: {t0/Ma2s} Myr")

  # plot the 1D velocity profile over time
  time_1d = np.linspace(0, 20, 201) * Ma2s # time array for plots
  bc_inv.plot_1D_velocity(time_1d)

  # material parameters
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

  # boundary conditions
  bc_phase_1 = boundary_conditions_Dirichlet(u)
  bc_phase_2 = boundary_conditions_Navier(u,grad_u,uL_2)
  
  # write the options for ptatin3d
  model_phase_1 = gp.Model(gp.InitialConditions(Domain,u),all_regions,bc_phase_1,output_path="extension")
  model_phase_2 = gp.Model(gp.InitialConditions(Domain,u),all_regions,bc_phase_2,output_path="compression")
  
  #with open("model_phase_1.sh","w") as f:
  #  f.write(model_phase_1.options)
  #with open("model_phase_2.sh","w") as f:
  #  f.write(model_phase_2.options)
  """
  time = np.linspace(0, 20, 21) * Ma2s # time array for plots
  root = "./"
  pvd = "timeseries.pvd"
  writer = gp.WriteVTS(Domain)
  bc_inv.paraview_velocity_inversion(writer,time,root,pvd)
  """

if __name__ == "__main__":
  main()