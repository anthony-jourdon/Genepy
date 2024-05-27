from bcpy import StokesBoundaryCondition

class Neumann(StokesBoundaryCondition):
  """
  .. py:class:: Neumann(tag:int, name:str, dev_stress_expression=None, mesh_file:str="path_to_file", model_name:str="model_GENE3D")

    Class to generate the options for `pTatin3d`_
    Neumann boundary condition for the Stokes problem.

    :param int tag: tag of the boundary condition given by gmsh
    :param str name: name of the boundary condition in the model (**optional**)
    :param str dev_stress_expression: expression for the deviatoric stress
    :param str mesh_file: path to the mesh file containing the facets of the boundary
    :param str model_name: name of the model (default: model_GENE3D)

    Attributes
    ----------

    .. py:attribute:: tag
      :type: int
      :canonical: bcpy.boundary_conditions.neumann.Neumann.tag

      Tag of the boundary condition

    .. py:attribute:: prefix
      :type: str
      :canonical: bcpy.boundary_conditions.neumann.Neumann.prefix

      Prefix: "bc_neumann" for the options

    .. py:attribute:: bc_type
      :type: int
      :canonical: bcpy.boundary_conditions.neumann.Neumann.bc_type

      Type of the boundary condition in the model. For Neumann: 1

    .. py:attribute:: bc_name
      :type: str
      :canonical: bcpy.boundary_conditions.neumann.Neumann.bc_name

      Name of the boundary condition in the model, name is arbitrary but providing one is mandatory

    .. py:attribute:: mesh_file
      :type: str
      :canonical: bcpy.boundary_conditions.neumann.Neumann.mesh_file

      Path to the mesh file containing the facets of the boundary

    .. py:attribute:: dev_stress
      :type: str
      :canonical: bcpy.boundary_conditions.neumann.Neumann.dev_stress

      Expression for the deviatoric stress. Can be None, in that case no deviatoric stress 
      is imposed, only the pressure.

    Methods
    -------
  """
  def __init__(self, tag:int, name:str, dev_stress_expression=None, mesh_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    StokesBoundaryCondition.__init__(self,tag,mesh_file,model_name)
    self.prefix     = "bc_neumann"
    self.bc_type    = 1
    self.bc_name    = name
    self.dev_stress = dev_stress_expression
  
  def sprint_option(self):
    """
    sprint_option(self)
    generates the string to be written in the options file for the Neumann boundary condition.
    
    :return: string to be written in the options file
    :rtype: str
    """
    s  = StokesBoundaryCondition.sprint_option(self)
    s += f"-{self.model_name}_poisson_pressure_apply\n"
    if self.dev_stress is not None:
      s += f"-{self.model_name}_{self.prefix}_dev_stress_{self.tag} {self.dev_stress} # deviatoric stress\n"
    return s
  
  def __str__(self):
    s = StokesBoundaryCondition.__str__(self)
    s += f"Pressure component is automatically calculated\n"
    s += f"\tDeviatoric stress: {self.dev_stress}\n"
    return s