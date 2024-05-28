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

    .. note:: 
      By default, when imposing a Neumann boundary condition such that the 
      imposed traction :math:`\\mathbf T` is given by
      
      .. math::
        \\mathbf T = \\boldsymbol \\sigma \\mathbf n
      
      with the **full** stress tensor

      .. math::
        \\boldsymbol \\sigma = \\boldsymbol \\tau - p \\boldsymbol I 

      it is necessary to provide *at least* a traction vector capable to 
      maintain the fluid (the rocks) inside the physical domain i.e. to 
      provide the **part of the pressure related to the density structure and
      the gravitational field**. For this purpose, `pTatin3d`_ computes

      .. math::
        \\nabla \\cdot \\nabla p_p = \\nabla \\cdot (\\rho \\mathbf g) 
      
      where :math:`\\rho` is the density, :math:`\\mathbf g` is the gravitational acceleration vector,
      and :math:`p_p` is refered as the Poisson pressure.

      If nothing more is provided, the Neumann boundary condition will be imposed as:

      .. math::
        \\mathbf T = -p_p \\mathbf n

      as :math:`\\mathbf n` points outward the domain.
      It is possible to add a scalar function (expression) :math:`\\tau(\\mathbf x,p_p,t)` 
      with the syntax ``x,y,z,p,t`` for the textual variable names 
      (any standard c math function can be used in the expression, see `this link <https://github.com/codeplea/tinyexpr>`_).
      In this case, the Neumann boundary condition will be imposed as:

      .. math::
        \\mathbf T = (\\tau - p_p) \\mathbf n


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