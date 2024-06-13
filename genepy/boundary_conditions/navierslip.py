#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: navierslip.py
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

class NavierSlip(StokesBoundaryCondition):
  """
  .. py:class:: NavierSlip(tag:int, name:str, grad_u, u_orientation, mesh_file:str="path_to_file", model_name:str="model_GENE3D")

    Class to generate the options for `pTatin3d`_
    using Navier-slip type of boundary conditions.

    :param int tag: tag of the boundary condition given by gmsh
    :param str name: name of the boundary condition in the model
    :param grad_u: velocity gradient to impose, expected shape: ``(3,3)``
    :param u_orientation: orientation of the velocity, expected shape: ``(1,3)``
    :param str mesh_file: path to the mesh file containing the facets of the boundary
    :param str model_name: name of the model (default: model_GENE3D)

    Attributes
    ----------

    .. py:attribute:: tag
      :type: int
      :canonical: genepy.boundary_conditions.navierslip.NavierSlip.tag

      Tag of the boundary condition

    .. py:attribute:: prefix 
      :type: str
      :canonical: genepy.boundary_conditions.navierslip.NavierSlip.prefix

      Prefix: "bc_navier" for the options

    .. py:attribute:: bc_type
      :type: int
      :canonical: genepy.boundary_conditions.navierslip.NavierSlip.bc_type

      Type of the boundary condition in the model. For Navier-slip: 6

    .. py:attribute:: bc_name
      :type: str
      :canonical: genepy.boundary_conditions.navierslip.NavierSlip.bc_name

      Name of the boundary condition in the model, name is arbitrary but providing one is mandatory

    .. py:attribute:: mesh_file
      :type: str
      :canonical: genepy.boundary_conditions.navierslip.NavierSlip.mesh_file

      Path to the mesh file containing the facets of the boundary

    Methods
    -------
  """
  def __init__(self, tag:int, name:str, grad_u, u_orientation, mesh_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,mesh_file,model_name)
    self.prefix     = "bc_navier"
    self.bc_type    = 6
    self.bc_name    = name
    self.grad_u     = grad_u
    self.uL         = u_orientation
    
  def sprint_option(self):
    """
    sprint_option(self)
    Returns the string to be added to the options file descibing the Navier-slip boundary condition.

    :return: string formatted for `pTatin3d`_ input file using `PETSc`_ options format
    :rtype: str
    """
    s = StokesBoundaryCondition.sprint_option(self)
    # collect the considered components of the gradient in a dictionnary
    gu = { "duxdx" : self.grad_u[0,0],
           "duxdz" : self.grad_u[0,2],
           "duzdx" : self.grad_u[2,0],
           "duzdz" : self.grad_u[2,2]}
    # remove spaces in the expressions
    gu = { k:self.format_expression(str(v)) for k,v in gu.items() }
    # check if the gradient is a constant or an expression
    if any(self.is_expression(gu[g]) for g in gu):
      s += f"-{self.model_name}_{self.prefix}_expression_{self.tag}\n"

    s += f"-{self.model_name}_{self.prefix}_penalty_{self.tag} 1.0e3 # penalty for Nitsche's method\n"
    for component,value in gu.items():
      s += f"-{self.model_name}_{self.prefix}_{component}_{self.tag} {value}\n"
    s += f"-{self.model_name}_{self.prefix}_uL_{self.tag} {self.uL[0]},{self.uL[2]} # velocity orientation on boundary\n"
    s += f"-{self.model_name}_{self.prefix}_mathcal_H_{self.tag} 0,1,0,1,1,1\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"\tVelocity gradient:\n"
    components = ['x','y','z']
    for i in range(3):
      for j in range(3):
        s += f"\t\tdu{components[i]}/d{components[j]}: {str(self.grad_u[i,j])}\n"
    s += f"\tOrientation of velocity: {self.uL[0]},{self.uL[2]}\n"
    return s
