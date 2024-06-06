#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: dirichlet.py
#
#  This file is part of Genepy.
#
#  Genepy is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  Genepy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with Genepy. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

from genepy import StokesBoundaryCondition

class Dirichlet(StokesBoundaryCondition):
  """
  .. py:class:: Dirichlet(tag:int, name:str, components, velocity, mesh_file:str="path_to_file", model_name:str="model_GENE3D")

    Class to generate the options for `pTatin3d`_ for a Dirichlet type boundary condition.
    Inherited from :py:class:`StokesBoundaryCondition`.

    :param int tag: tag of the boundary condition given by gmsh
    :param str name: name of the boundary condition in the model
    :param list components: list of string of the velocity components on which to impose the Dirichlet condition. Possible values ``["x","y","z"]``
    :param velocity: velocity to impose, expected shape: ``(1,dim)``
    :param str mesh_file: path to the mesh file containing the facets of the boundary
    :param str model_name: name of the model (default: model_GENE3D)

    Attributes
    ----------

    .. py:attribute:: tag
      :type: int
      :canonical: genepy.boundary_conditions.dirichlet.Dirichlet.tag

      Tag of the boundary condition
    
    .. py:attribute:: prefix
      :type: str
      :canonical: genepy.boundary_conditions.dirichlet.Dirichlet.prefix

      Prefix: "bc_dirichlet" for the options

    .. py:attribute:: bc_type
      :type: int
      :canonical: genepy.boundary_conditions.dirichlet.Dirichlet.bc_type

      Type of the boundary condition in the model. For Dirichlet: 7
    
    .. py:attribute:: bc_name
      :type: str
      :canonical: genepy.boundary_conditions.dirichlet.Dirichlet.bc_name

      Name of the boundary condition in the model, name is arbitrary but providing one is mandatory
    
    .. py:attribute:: mesh_file
      :type: str
      :canonical: genepy.boundary_conditions.dirichlet.Dirichlet.mesh_file

      Path to the mesh file containing the facets of the boundary
  """
  def __init__(self, tag:int, name:str, components, velocity, mesh_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,mesh_file,model_name)
    self.prefix     = "bc_dirichlet"
    self.bc_type    = 7
    self.bc_name    = name
    self.u          = velocity
    self.components = components
  
  def sprint_option(self):
    """
    sprint_option(self)
    Returns the string to be added to the options file descibing the Dirichlet boundary condition.
    Calls :meth:`genepy.boundary_conditions.bcs.StokesBoundaryCondition.sprint_option` first.

    :return: string to be added to the options file
    :rtype: str
    """
    s = StokesBoundaryCondition.sprint_option(self)
    for d in self.components:
      if   d == "x": dim = 0
      elif d == "y": dim = 1
      elif d == "z": dim = 2
      u_string = str(self.u[0,dim])
      u_split  = u_string.split()
      nmembers = len(u_split)
      u_nospace = ""
      for j in range(nmembers):
        u_nospace += u_split[j]
      s += f"-{self.model_name}_{self.prefix}_u{d}_{self.tag} {u_nospace}\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"\tComponents: {self.components}\n"
    s += f"\tVelocity:\n"
    for d in self.components:
      if   d == "x": dim = 0
      elif d == "y": dim = 1
      elif d == "z": dim = 2
      s += f"\t\tu{d}: {self.u[0,dim]}\n"
    return s

class DirichletUdotN(StokesBoundaryCondition):
  """
  .. py:class:: DirichletUdotN(tag:int, name:str, mesh_file:str="path_to_file", model_name:str="model_GENE3D")

    Class to generate the options for `pTatin3d`_ for a special Dirichlet type boundary condition.
    This boundary condition is used to impose a velocity on the basal face of the mesh such that the inflow and outflow of the entire domain are balanced.

    :param int tag: tag of the boundary condition given by gmsh
    :param str name: name of the boundary condition in the model
    :param str mesh_file: path to the mesh file containing the facets of the boundary
    :param str model_name: name of the model (default: model_GENE3D)

    Attributes
    ----------

    .. py:attribute:: tag
      :type: int
      :canonical: genepy.boundary_conditions.dirichlet.DirichletUdotN.tag

      Tag of the boundary condition

    .. py:attribute:: prefix
      :type: str
      :canonical: genepy.boundary_conditions.dirichlet.DirichletUdotN.prefix

      Prefix: "bc_dirichlet" for the options

    .. py:attribute:: bc_type
      :type: int
      :canonical: genepy.boundary_conditions.dirichlet.DirichletUdotN.bc_type

      Type of the boundary condition in the model. For Dirichlet: 7

    .. py:attribute:: bc_name
      :type: str
      :canonical: genepy.boundary_conditions.dirichlet.DirichletUdotN.bc_name

      Name of the boundary condition in the model, name is arbitrary but providing one is mandatory

    .. py:attribute:: mesh_file
      :type: str
      :canonical: genepy.boundary_conditions.dirichlet.DirichletUdotN.mesh_file

      Path to the mesh file containing the facets of the boundary

    Methods
    -------
  """
  def __init__(self, tag:int, name:str, mesh_file:str="path_to_file",model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,mesh_file,model_name)
    self.prefix     = "bc_dirichlet"
    self.bc_type    = 7
    self.bc_name    = name
    
  def sprint_option(self):
    """
    sprint_option(self)
    Returns the string to be added to the options file descibing the special u.n Dirichlet boundary condition.
    Calls :meth:`genepy.boundary_conditions.bcs.StokesBoundaryCondition.sprint_option` first.

    :return: string to be added to the options file
    :rtype: str
    """
    s = StokesBoundaryCondition.sprint_option(self)
    s += f"-{self.model_name}_{self.prefix}_bot_u.n_{self.tag}\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"\tVelocity: balancing inflow/outflow\n"
    s += f"\ton bottom face using int u.n dS\n"
    return s
