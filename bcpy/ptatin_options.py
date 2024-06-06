#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: ptatin_options.py
#
#  This file is part of bc-pre-processing.
#
#  bc-pre-processing is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  bc-pre-processing is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with bc-pre-processing. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

from bcpy import ModelRegions
from bcpy import ModelBCs
from bcpy import Domain
from bcpy import MarkersManagement
from bcpy import InitialConditions
import numpy as np

class Model:
  """
  .. py:class:: Model(model_ics:InitialConditions, model_regions:ModelRegions, model_bcs:ModelBCs, model_name:str="model_GENE3D", **kwargs)

    Class to define the options for `pTatin3d`_ model.

    :param InitialConditions model_ics: Initial conditions of the model.
    :param ModelRegions model_regions: Material parameters of the model.
    :param ModelBCs model_bcs: Boundary conditions of the model.
    :param str model_name: Name of the model. Default is ``"model_GENE3D"``.
    :param kwargs: Additional arguments for the model.

    **Keyword arguments (optional)**

    :param int output_frequency: Frequency of results output. Default is ``25``.
    :param str output_path: Path of the output directory. Default is ``"output"``.
    :param int checkpoint_ncpumins: Checkpoint every n cpu minutes. Default is ``230``.
    :param int checkpoint_frequency: Checkpoint every n time step. Default is ``25``.
    :param float length_bar: Length scaling. Default is ``1.0e5`` m.
    :param float viscosity_bar: Viscosity scaling. Default is :math:`10^{22}` Pa.s.
    :param float velocity_bar: Velocity scaling. Default is :math:`10^{-10}` m.s\ :sup:`-1`.
    :param int nsteps: Max number of time steps. Default is ``100000``.
    :param float dtmin: Min time step. Default is ``1.0e-6``.
    :param float dtmax: Max time step. Default is ``0.5``.
    :param float isostatic_density: Density reference for initial isostatic topography. Default is ``3300`` kg.m\ :sup:`-3`.
    :param float isostatic_depth: Depth of compensation for initial isostatic topography. Default is ``-40`` km.
    :param float max_surface_displacement: Max surface displacement per time step. Default is 500 m.
    :param MarkersManagement markers: Markers layout and population control options.
    :param float viscosity_min: Minimum viscosity value. Default is :math:`10^{19}` Pa.s.
    :param float viscosity_max: Maximum viscosity value. Default is :math:`10^{25}` Pa.s.
    :param SPMDiffusion spm: Surface processes SPMDiffusion class instance.
    :param Pswarm pswarm: Passive tracers class instance.
    :param int mpi_ranks: Number of MPI ranks, it is used to define the solver options. 
                          However, only ``1`` and ``1024`` MPI ranks are actually supported.
                          Any value other than ``1`` will set the solver options for ``1024`` MPI ranks.
                          Default is ``1``. The solver options for ``1024`` MPI ranks 
                          may not be the best for a particular problem.
                          
    Attributes
    ----------

    .. py:attribute:: name
      :type: str

        Name of the model. Default is ``"model_GENE3D"``.

    .. py:attribute:: regions
      :type: ModelRegions

        Material parameters of the model.

    .. py:attribute:: bcs
      :type: ModelBCs

        Boundary conditions of the model.

    .. py:attribute:: ics
      :type: InitialConditions

        Initial conditions of the model.

    .. py:attribute:: options
      :type: str

        Options for the model.

    Example

  """
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

    # initial topography options
    density = kwargs.get("isostatic_density", 3300.0)
    depth   = kwargs.get("isostatic_depth", -40.0e3)
    self.isostatic_topography(density, depth)
    
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

    # time stepping options
    nsteps = kwargs.get("nsteps", 100000)
    dtmin  = kwargs.get("dtmin", 1.0e-6)
    dtmax  = kwargs.get("dtmax", 0.5)
    max_surface_displacement = kwargs.get("max_surface_displacement", 500.0) / length_bar
    self.time_stepping(nsteps, dtmin, dtmax, max_surface_displacement)
    
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
  
  def time_stepping(self, nsteps:int, dtmin:float, dtmax:float, max_surface_displacement:float) -> None:
    self.options += "########### Time Stepping ###########\n"
    self.options += f"-nsteps {nsteps} # max number of time steps\n"
    self.options += f"-dt_min {dtmin} # min time step\n"
    self.options += f"-dt_max {dtmax} # max time step\n"
    self.options += f"-dt_max_surface_displacement {max_surface_displacement} # max surface displacement per time step\n"
    return
  
  def isostatic_topography(self, density:float, depth:float):
    prefix = "isostatic"
    self.options += f"########### {prefix} initial topography ###########\n"
    self.options += f"-{self.name}_{prefix}_density_ref {density}\n"
    self.options += f"-{self.name}_{prefix}_depth {depth}\n"
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
