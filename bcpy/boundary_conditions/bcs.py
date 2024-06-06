#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: bcs.py
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

class BoundaryCondition:
  def __init__(self, model_name:str="model_GENE3D") -> None:
    self.model_name = model_name

class StokesBoundaryCondition(BoundaryCondition):
  """
  .. py:class:: StokesBoundaryCondition(tag:int=0, mesh_file:str="path_to_file", model_name:str="model_GENE3D")

    Class to generate the options for `pTatin3d`_ 
    boundary condition for the Stokes problem.

    :param int tag: tag of the boundary condition given by gmsh
    :param str mesh_file: path to the mesh file containing the facets of the boundary
    :param str model_name: name of the model (default: model_GENE3D)

    Attributes
    ----------

    .. py:attribute:: tag
      :type: int
      :canonical: bcpy.boundary_conditions.bcs.StokesBoundaryCondition.tag

      Tag of the boundary condition

    .. py:attribute:: prefix
      :type: str
      :canonical: bcpy.boundary_conditions.bcs.StokesBoundaryCondition.prefix

      Prefix: "bc" for the options

    .. py:attribute:: bc_type
      :type: int
      :canonical: bcpy.boundary_conditions.bcs.StokesBoundaryCondition.bc_type

      Type of the boundary condition in the model

    .. py:attribute:: bc_name
      :type: str
      :canonical: bcpy.boundary_conditions.bcs.StokesBoundaryCondition.bc_name

      Name of the boundary condition in the model

    .. py:attribute:: mesh_file
      :type: str
      :canonical: bcpy.boundary_conditions.bcs.StokesBoundaryCondition.mesh_file

      Path to the mesh file containing the facets of the boundary
  """
  def __init__(self, tag:int=0, mesh_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    BoundaryCondition.__init__(self,model_name)
    self.tag:int        = tag
    self.prefix:str     = "bc"
    self.bc_type:int    = 0
    self.bc_name:str    = "uninitialized"
    self.mesh_file:str  = mesh_file
  
  def sprint_option(self):
    """
    sprint_option(self)
    Returns a string formatted for `pTatin3d`_ input file using `PETSc`_ options format.
    All subclasses of :class:`StokesBoundaryCondition` call this method.

    :return: string formatted for `pTatin3d`_ input file
    :rtype: str
    """
    s  = f"###### Boundary condition {self.bc_name} {self.tag} ######\n"
    s += f"###### {self.__class__.__name__} boundary conditions ######\n"
    s += f"-{self.model_name}_bc_sc_name_{self.tag} {self.bc_name} # name of the boundary (arbitrary)\n"
    s += f"-{self.model_name}_bc_sc_type_{self.tag} {self.bc_type} # type of BC\n"
    s += f"-{self.model_name}_bc_facet_mesh_file_{self.tag} {self.mesh_file} # path to file of the boundary mesh\n"
    return s

  def __str__(self):
    s = f"{self.__class__.__name__}:\n"
    s += f"\tSurface tag: {self.tag}\n"
    s += f"\tName:        {self.bc_name}\n"
    s += f"\tType:        {self.bc_type}\n"
    return s

class TemperatureBC(BoundaryCondition):
  def __init__(self, conditions:dict, model_name:str="model_GENE3D") -> None:
    BoundaryCondition.__init__(self,model_name)
    self.conditions = conditions
  
  def sprint_option(self):
    s = f"###### Temperature boundary conditions ######\n"
    for face in self.conditions:
      if face not in ['xmin','xmax','ymin','ymax','zmin','zmax']:
        err = f"Error: face {face} not recognized.\n"
        err += "Recognized faces are: xmin, xmax, ymin, ymax, zmin, zmax.\n"
        raise ValueError(err)
      s += f"-{self.model_name}_bc_energy_{face} {self.conditions[face]} # Temperature BC on face {face}\n"
    return s
  
  def __str__(self):
    s = "Temperature boundary conditions:\n"
    for face,value in zip(self.faces,self.values):
      s += f"\tFace: {face}\n"
      s += f"\tValue: {value}\n"
    return s

class ModelBCs:
  """
  .. py:class:: ModelBCs(velocity_bcs:list[StokesBoundaryCondition], energy_bc:TemperatureBC=None, model_name:str="model_GENE3D")

    Class to generate options for the boundary conditions of a `pTatin3d`_ model.
    Enforce the model name for all boundary conditions objects.

    :param list[StokesBoundaryCondition] velocity_bcs: list of instances of class 
      :py:class:`StokesBoundaryCondition` 
      for the Stokes boundary conditions of the model.
    :param TemperatureBC energy_bc: instance of class :py:class:`TemperatureBC` 
      for the energy boundary conditions of the model.
    :param str model_name: name of the model (default: model_GENE3D)

    Attributes
    ----------
    .. py:attribute:: model_name
      :type: str
      :canonical: bcpy.boundary_conditions.bcs.ModelBCs.model_name

      Name of the model

    .. py:attribute:: u_bcs
      :type: list[StokesBoundaryCondition]
      :canonical: bcpy.boundary_conditions.bcs.ModelBCs.u_bcs

      List of instances of class :py:class:`StokesBoundaryCondition` (and its children) 
      for the Stokes problem boundary conditions.

    .. py:attribute:: energy_bc
      :type: TemperatureBC
      :canonical: bcpy.boundary_conditions.bcs.ModelBCs.energy_bc

      Instance of class :py:class:`TemperatureBC` for the energy boundary conditions

    Example
    -------

    The following example assumes that the velocity function ``u``, its gradient ``grad_u``
    and the orientation vector ``uL`` are defined elsewhere in the code.

    .. code-block:: python

      import bcpy as bp
    
      # Velocity boundary conditions
      bcs = [
          bp.Dirichlet(23,"Zmax",["x","z"],u),
          bp.Dirichlet(37,"Zmin",["x","z"],u),
          bp.NavierSlip(32,"Xmax",grad_u,uL),
          bp.NavierSlip(14,"Xmin",grad_u,uL),
          bp.DirichletUdotN(33,"Bottom"),
      ]
      # Temperature boundary conditions
      Tbcs = bp.TemperatureBC({"ymax":0.0, "ymin":1450.0})
      # collect all boundary conditions
      all_bcs = bp.ModelBCs(bcs,Tbcs)
  """
  def __init__(self, velocity_bcs:list[StokesBoundaryCondition], 
               energy_bc:TemperatureBC=None, 
               model_name:str="model_GENE3D") -> None:
    self.model_name = model_name
    self.u_bcs      = velocity_bcs
    self.energy_bc  = energy_bc
    # enforce model name for all boundary conditions
    for bc in self.u_bcs:
      bc.model_name = self.model_name
    if self.energy_bc is not None:
      self.energy_bc.model_name = self.model_name

  def add_velocity_bc(self, bc:StokesBoundaryCondition) -> None:
    self.u_bcs.append(bc)
  
  def sprint_option(self):
    s = "########### Boundary conditions ###########\n"
    nsubfaces = len(self.u_bcs)
    s += f"-{self.model_name}_bc_nsubfaces {nsubfaces} # number of subfaces\n"
    s += f"-{self.model_name}_bc_tag_list "
    for n in range(nsubfaces-1):
      s += f"{self.u_bcs[n].tag},"
    s += f"{self.u_bcs[nsubfaces-1].tag} # list of boundary tags, determine the order of application\n"
    for bc in self.u_bcs:
      s += bc.sprint_option()
    if self.energy_bc is not None:
      s += self.energy_bc.sprint_option()
    return s

  def __str__(self):
    s = f"{self.model_name} boundary conditions:\n"
    for bc in self.u_bcs:
      s += str(bc)
    if self.energy_bc is not None:
      s += str(self.energy_bc)
    return s
