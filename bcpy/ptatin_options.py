from bcpy import ModelRegions
from bcpy import ModelBCs
from bcpy import Domain
from bcpy import MarkersManagement
from bcpy import InitialConditions
import numpy as np

class Model:
  def __init__(self, model_ics:InitialConditions, 
               model_regions:ModelRegions, 
               model_bcs:ModelBCs,
               model_name:str="model_GENE3D",
               **kwargs) -> None:
    self.name    = model_name
    self.regions = model_regions
    self.bcs     = model_bcs
    self.ics     = model_ics

    self.options  = "########### Gene3D options file ###########\n"
    self.options += "-ptatin_model Gene3D\n"
    self.options += "-activate_energyfv true # use finite volume for temperature (Do not change)\n"
    self.options += "-view_ic # write .vts file for initial temperature\n"
    self.options += f"-{self.name}_ic_temperature_from_file # load temperature from petsc vec file\n"
    
    # output options
    output_frequency = kwargs.get("output_frequency", 25)
    output_path      = kwargs.get("output_path", "output")
    self.output(output_frequency, output_path)

    # checkpointing options
    ncpumins   = kwargs.get("checkpoint_ncpumins", 230)
    chkpt_freq = kwargs.get("checkpoint_frequency", output_frequency)
    self.checkpointing(ncpumins, chkpt_freq)

    # initial conditions options
    self.ics.model_name = self.name
    self.options += self.ics.sprint_option()
    
    # spm options
    if "spm" in kwargs:
      kwargs["spm"].model_name = self.name
      self.options += kwargs["spm"].sprint_option()
    
    # markers layout and population control options
    markers = kwargs.get("markers", MarkersManagement())
    self.options += markers.sprint_option()

    # scaling options
    length_bar    = kwargs.get("length_bar",    1.0e5)
    viscosity_bar = kwargs.get("viscosity_bar", 1.0e22)
    velocity_bar  = kwargs.get("velocity_bar",  1.0e-10)
    self.scaling(length_bar, viscosity_bar, velocity_bar)
    
    # passive swarm options
    if "pswarm" in kwargs:
      pswarm = kwargs["pswarm"]
      pswarm.model_name = self.name
      attributes = vars(pswarm)
      # scale length because pswarm expect scaled values
      if "O" in attributes:
        pswarm.O = np.asarray(pswarm.O, dtype=np.float32) / length_bar
      if "L" in attributes:
        pswarm.L = np.asarray(pswarm.L, dtype=np.float32) / length_bar
      self.options += pswarm.sprint_option()

    # boundary conditions options
    self.bcs.model_name = self.name
    self.options += self.bcs.sprint_option()

    # viscosity cutoff options
    viscosity_min = kwargs.get("viscosity_min", 1.e+19)
    viscosity_max = kwargs.get("viscosity_max", 1.e+25)
    self.viscosity_cutoff(viscosity_min, viscosity_max)

    # material parameters options
    self.regions.model_name = self.name
    self.options += self.regions.sprint_option()

    # solvers options
    ranks = kwargs.get("mpi_ranks", 1)
    if ranks > 1: self.solvers_1024MPIranks()
    else:         self.solvers_1MPIrank()

    return

  def output(self, output_frequency:int, output_path:str) -> None:
    prefix = "output"
    self.options +=f"########### {prefix} ###########\n"
    self.options +=f"-{prefix}_frequency {output_frequency} # frequency of results output\n"
    self.options +=f"-{prefix}_path {output_path} # path of the output directory\n"
    return
  
  def viscosity_cutoff(self, viscosity_min:float, viscosity_max:float) -> None:
    if viscosity_min > viscosity_max:
      err = f"Error: viscosity_min ({viscosity_min}) > viscosity_max ({viscosity_max}).\n"
      raise ValueError(err)
    
    prefix = "viscosity_cutoff"
    self.options +=f"########### {prefix} ###########\n"
    self.options +=f"-{self.name}_{prefix}_apply\n"
    self.options +=f"-{self.name}_{prefix}_eta_lower {viscosity_min}\n"
    self.options +=f"-{self.name}_{prefix}_eta_upper {viscosity_max}\n"
    return 
  
  def checkpointing(self, ncpumins:int, chkpt_freq:int) -> None:
    prefix = "checkpoint"
    self.options +=f"########### {prefix} ###########\n"
    self.options +=f"-{prefix}_every_ncpumins {ncpumins} # checkpoint every n cpu minutes\n"
    self.options +=f"-{prefix}_every_nsteps {chkpt_freq} # checkpoint every n time step\n"
    return

  def scaling(self,length_bar:float, viscosity_bar:float, velocity_bar:float) -> None:
    prefix = "scaling"
    self.options += f"########### {prefix} ###########\n"
    self.options += f"-{self.name}_{prefix}_length {length_bar} # length scaling\n"
    self.options += f"-{self.name}_{prefix}_viscosity {viscosity_bar} # viscosity scaling\n"
    self.options += f"-{self.name}_{prefix}_velocity {velocity_bar} # velocity scaling\n"
    return
  
  def solvers_1MPIrank(self) -> None:
    self.options += f"-{self.name}_output_markers\n"
    self.options += f"-{self.name}_bc_debug\n"
    self.options += "##########################################\n"
    self.options += "########### SOLVERS 1 MPI rank ###########\n"
    self.options += "##########################################\n"
    
    self.options += "########### Stokes ###########\n"
    self.options += "-A11_operator_type 2,2,0\n"
    self.options += "-dau_nlevels 3\n"
    self.options += "-stk_velocity_dm_mat_type aij\n"
    
    self.options += "###### SNES ######\n"
    self.options += "-newton_its 0\n"
    self.options += "-snes_atol 1e-6\n"
    self.options += "-snes_max_it 5\n"
    self.options += "-snes_rtol 1e-8\n"
    self.options += "-snes_monitor\n"
    self.options += "-snes_converged_reason\n"
    self.options += "-snes_linesearch_type basic\n"
    
    self.options += "###### KSP ######\n"
    self.options += "-ksp_max_it 200\n"
    self.options += "-ksp_rtol 1e-3\n"
    self.options += "-ksp_type fgmres\n"
    self.options += "-pc_type fieldsplit\n"
    self.options += "-pc_fieldsplit_schur_fact_type upper\n"
    self.options += "-pc_fieldsplit_type schur\n"
    
    self.options += "\t###### Fieldsplit p ######\n"
    self.options += "\t-fieldsplit_p_ksp_type preonly\n"
    self.options += "\t-fieldsplit_p_pc_type bjacobi\n"
    
    self.options += "\t###### Fieldsplit u ######\n"
    self.options += "\t-fieldsplit_u_ksp_max_it 200\n"
    self.options += "\t-fieldsplit_u_ksp_monitor\n"
    self.options += "\t-fieldsplit_u_ksp_converged_reason\n"
    self.options += "\t-fieldsplit_u_ksp_rtol 1.0e-2\n"
    self.options += "\t-fieldsplit_u_ksp_type fgmres\n"
    self.options += "\t-fieldsplit_u_pc_type mg\n"
    self.options += "\t-fieldsplit_u_pc_mg_levels 3\n"
    self.options += "\t-fieldsplit_u_pc_mg_cycle_type v\n"
    self.options += "\t-fieldsplit_u_mg_levels_ksp_type fgmres\n"
    self.options += "\t-fieldsplit_u_mg_levels_pc_type bjacobi\n"
    self.options += "\t-fieldsplit_u_mg_levels_ksp_max_it 8\n"
    self.options += "\t-fieldsplit_u_mg_levels_ksp_norm_type NONE\n"
    
    self.options += "\t\t###### Coarse Grid ######\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_pc_type cholesky\n"
    
    
    self.options += "########### Temperature ###########\n"
    self.options += "-ptatin_energyfv_nsub 3,3,3\n"
    self.options += "-energyfv_snes_converged_reason\n"
    self.options += "-energyfv_ksp_type fgmres\n"
    self.options += "-energyfv_snes_ksp_ew\n"
    self.options += "-energyfv_snes_ksp_ew_rtol0 1.0e-2\n"
    self.options += "-energyfv_snes_ksp_ew_version 1\n"
    self.options += "-energyfv_snes_lag_preconditioner -2\n"
    self.options += "-energyfv_snes_monitor\n"
    self.options += "-energyfv_snes_rtol 1.0e-6\n"
    self.options += "-energyfv_mg_levels_esteig_ksp_norm_type none\n"
    self.options += "-energyfv_mg_levels_esteig_ksp_type cg\n"
    self.options += "-energyfv_mg_levels_ksp_chebyshev_esteig 0,0.01,0,1.1\n"
    self.options += "-energyfv_mg_levels_ksp_max_it 8\n"
    self.options += "-energyfv_mg_levels_ksp_type gmres\n"
    self.options += "-energyfv_mg_levels_pc_type bjacobi\n"
    self.options += "-energyfv_pc_type gamg\n"
    self.options += "-fvpp_ksp_monitor\n"
    self.options += "-fvpp_ksp_rtol 1.0e-10\n"
    self.options += "-fvpp_ksp_type cg\n"
    self.options += "-fvpp_mg_coarse_pc_type redundant\n"
    self.options += "-fvpp_mg_coarse_redundant_pc_factor_mat_solver_type gamg\n"
    self.options += "-fvpp_mg_coarse_redundant_pc_type gamg\n"
    self.options += "-fvpp_mg_levels_esteig_ksp_norm_type none\n"
    self.options += "-fvpp_mg_levels_esteig_ksp_type cg\n"
    self.options += "-fvpp_mg_levels_ksp_chebyshev_esteig 0,0.01,0,1.1\n"
    self.options += "-fvpp_mg_levels_ksp_max_it 4\n"
    self.options += "-fvpp_mg_levels_ksp_norm_type none\n"
    self.options += "-fvpp_mg_levels_ksp_type chebyshev\n"
    self.options += "-fvpp_mg_levels_pc_type jacobi\n"
    self.options += "-fvpp_operator_fvspace false\n"
    self.options += "-fvpp_pc_type gamg\n"


    self.options += "########### Poisson Pressure ###########\n"
    self.options += "-LP_snes_rtol 1.0e-6\n"
    self.options += "-LP_snes_monitor\n"
    self.options += "-LP_snes_converged_reason\n"
    self.options += "-LP_ksp_monitor\n"
    self.options += "-LP_ksp_converged_reason\n"
    self.options += "-LP_ksp_atol 1.0e-10\n"
    self.options += "-LP_ksp_rtol 1.0e-6\n"
    self.options += "-LP_ksp_type fgmres\n"
    self.options += "-LP_pc_type bjacobi\n"
    self.options += "-LP_snes_ksp_ew\n"
    self.options += "-LP_snes_ksp_ew_rtol0 1.0e-2\n"
    self.options += "-LP_snes_ksp_ew_version 1\n"
    self.options += "-LP_snes_lag_preconditioner -2\n"
    return
  
  def solvers_1024MPIranks(self) -> None:
    self.options += "##############################################\n"
    self.options += "########### SOLVERS 1024 MPI ranks ###########\n"
    self.options += "##############################################\n"

    self.options += "########### Stokes ###########\n"
    self.options += "-A11_operator_type 2,2,0,1\n"
    self.options += "-dau_nlevels 4\n"
    self.options += "-stk_velocity_dm_mat_type aij\n"
    self.options += "-a11_op avx\n"

    self.options += "###### SNES ######\n"
    self.options += "-newton_its 0\n"
    self.options += "-snes_atol 1e-6\n"
    self.options += "-snes_max_it 5\n"
    self.options += "-snes_rtol 1e-3\n"
    self.options += "-snes_monitor\n"
    self.options += "-snes_converged_reason\n"
    self.options += "-snes_linesearch_type basic\n"

    self.options += "###### KSP ######\n"
    self.options += "-ksp_max_it 200\n"
    self.options += "-ksp_rtol 1e-3\n"
    self.options += "-ksp_type fgmres\n"
    self.options += "-pc_type fieldsplit\n"
    self.options += "-pc_fieldsplit_schur_fact_type upper\n"
    self.options += "-pc_fieldsplit_type schur\n"

    self.options += "\t###### Fieldsplit p ######\n"
    self.options += "\t-fieldsplit_p_ksp_type preonly\n"
    self.options += "\t-fieldsplit_p_pc_type jacobi\n"

    self.options += "\t###### Fieldsplit u ######\n"
    self.options += "\t-fieldsplit_u_ksp_max_it 200\n"
    self.options += "\t-fieldsplit_u_ksp_monitor\n"
    self.options += "\t-fieldsplit_u_ksp_converged_reason\n"
    self.options += "\t-fieldsplit_u_ksp_rtol 1.0e-4\n"
    self.options += "\t-fieldsplit_u_ksp_type fgmres\n"
    self.options += "\t-fieldsplit_u_pc_type mg\n"
    self.options += "\t-fieldsplit_u_pc_mg_levels 4\n"
    self.options += "\t-fieldsplit_u_mg_levels_ksp_type fgmres\n"
    self.options += "\t-fieldsplit_u_mg_levels_pc_type jacobi\n"
    self.options += "\t-fieldsplit_u_mg_levels_ksp_max_it 10\n"
    self.options += "\t-fieldsplit_u_mg_levels_ksp_norm_type NONE\n"

    self.options += "\t\t###### Coarse Grid ######\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_pc_type telescope\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_ksp_type preonly\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_pc_telescope_reduction_factor 16\n"
    self.options += "\t\t# From 1024 MPI ranks, new there are 1024/16 = 64 MPI ranks\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_pc_type mg\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_ksp_type preonly\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_pc_mg_levels 3\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_pc_mg_galerkin\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_mg_levels_ksp_type fgmres\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_mg_levels_pc_type jacobi\n"
    self.options += "\t\t-fieldsplit_u_mg_coarse_telescope_mg_levels_ksp_max_it 10\n"
    
    self.options += "\t\t\t###### Coarse Coarse Grid ######\n"
    self.options += "\t\t\t-fieldsplit_u_mg_coarse_telescope_mg_coarse_pc_type telescope\n"
    self.options += "\t\t\t-fieldsplit_u_mg_coarse_telescope_mg_coarse_pc_telescope_reduction_factor 64\n"
    self.options += "\t\t\t-fieldsplit_u_mg_coarse_telescope_mg_coarse_telescope_pc_type cholesky\n"

    self.options += "########### Temperature ###########\n"
    self.options += "-ptatin_energyfv_nsub 1,1,1\n"
    self.options += "-energyfv_snes_converged_reason\n"
    self.options += "-energyfv_snes_ksp_ew\n"
    self.options += "-energyfv_snes_ksp_ew_rtol0 1.0e-2\n"
    self.options += "-energyfv_snes_ksp_ew_version 1\n"
    self.options += "-energyfv_snes_lag_preconditioner -2\n"
    self.options += "-energyfv_snes_monitor\n"
    self.options += "-energyfv_snes_rtol 1.0e-6\n"
    self.options += "-energyfv_snes_atol 1.0e-6\n"
    self.options += "-energyfv_ksp_monitor\n"
    self.options += "-energyfv_ksp_type fgmres\n"
    self.options += "-energyfv_pc_type gamg\n"
    self.options += "\t-energyfv_mg_levels_ksp_max_it 8\n"
    self.options += "\t-energyfv_mg_levels_ksp_type fgmres\n"
    self.options += "\t-energyfv_mg_levels_pc_type jacobi\n"
    self.options += "\t-energyfv_mg_coarse_pc_type telescope\n"
    self.options += "\t-energyfv_mg_coarse_pc_telescope_reduction_factor 16\n"
    self.options += "\t# MPI ranks: 1024/16 = 64\n"
    self.options += "\t\t-energyfv_mg_coarse_telescope_pc_type gamg\n"
    self.options += "\t\t-energyfv_mg_coarse_telescope_mg_levels_ksp_type fgmres\n"
    self.options += "\t\t-energyfv_mg_coarse_telescope_mg_levels_pc_type jacobi\n"
    self.options += "\t\t-energyfv_mg_coarse_telescope_mg_levels_ksp_max_it 8\n"
    self.options += "\t\t\t# Coarse Coarse\n"
    self.options += "\t\t\t-energyfv_mg_coarse_telescope_mg_coarse_pc_type telescope\n"
    self.options += "\t\t\t-energyfv_mg_coarse_telescope_mg_coarse_pc_telescope_reduction_factor 64\n"
    self.options += "\t\t\t-energyfv_mg_coarse_telescope_mg_coarse_telescope_pc_type lu\n"
    self.options += "-fvpp_ksp_monitor\n"
    self.options += "-fvpp_ksp_rtol 1.0e-10\n"
    self.options += "-fvpp_ksp_type cg\n"
    self.options += "-fvpp_mg_coarse_pc_type redundant\n"
    self.options += "-fvpp_mg_coarse_redundant_pc_factor_mat_solver_type gamg\n"
    self.options += "-fvpp_mg_coarse_redundant_pc_type gamg\n"
    self.options += "-fvpp_mg_levels_esteig_ksp_norm_type none\n"
    self.options += "-fvpp_mg_levels_esteig_ksp_type cg\n"
    self.options += "-fvpp_mg_levels_ksp_chebyshev_esteig 0,0.01,0,1.1\n"
    self.options += "-fvpp_mg_levels_ksp_max_it 4\n"
    self.options += "-fvpp_mg_levels_ksp_norm_type none\n"
    self.options += "-fvpp_mg_levels_ksp_type chebyshev\n"
    self.options += "-fvpp_mg_levels_pc_type jacobi\n"
    self.options += "-fvpp_operator_fvspace false\n"
    self.options += "-fvpp_pc_type gamg\n"

    self.options += "########### Poisson Pressure ###########\n"
    self.options += "-LP_snes_rtol 1.0e-6\n"
    self.options += "-LP_snes_monitor\n"
    self.options += "-LP_snes_converged_reason\n"
    self.options += "-LP_ksp_monitor\n"
    self.options += "-LP_ksp_converged_reason\n"
    self.options += "-LP_ksp_atol 1.0e-10\n"
    self.options += "-LP_ksp_rtol 1.0e-6\n"
    self.options += "-LP_ksp_type fgmres\n"
    self.options += "-LP_pc_type bjacobi\n"
    self.options += "-LP_snes_ksp_ew\n"
    self.options += "-LP_snes_ksp_ew_rtol0 1.0e-2\n"
    self.options += "-LP_snes_ksp_ew_version 1\n"
    self.options += "-LP_snes_lag_preconditioner -2\n"
    return