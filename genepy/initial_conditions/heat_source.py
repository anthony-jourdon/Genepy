from genepy.initial_conditions.gaussian import Gaussian

class InitialHeatSource:
  """
  .. py:class:: InitialHeatSource(hs_expression)

    Class to generate options for the initial heat source when using the class
    :py:class:`genepy.EnergySourceMaterialPointValue`. 
    The heat source should be given as an expression or an instance of the 
    :py:class:`genepy.Gaussian` class.
    If the expression is a constant value, prefer to use the
    :py:class:`genepy.EnergySourceConstant` class.

    Example:
    --------
    In the following example we assume that the :py:class:`Domain <genepy.Domain>` class is already instantiated.

    .. code-block:: python

      import numpy as np
      import genepy as gp

      # shape of the gaussian
      coeff = 0.5 * 6.0e-5**2
      a = np.array([coeff],dtype=np.float64)
      b = np.zeros(shape=(1),dtype=np.float64)
      c = np.array([coeff],dtype=np.float64)

      # position of the centre of the gaussian
      x0 = np.array([0.5 * (Domain.L_num[0] - Domain.O_num[0])], dtype=np.float64)
      z0 = np.array([0.5 * (Domain.L_num[2] - Domain.O_num[2])], dtype=np.float64)

      # amplitude of the gaussian (heat source)
      A = np.array([1.5e-6],dtype=np.float64)

      # Create gaussian object
      Gaussian = gp.Gaussian(Domain,1,A,a,b,c,x0,z0)
      Gaussian.evaluate_gaussians()

      # Create initial heat source
      Hini = gp.InitialHeatSource(Gaussian)

    :param hs_expression: heat source expression (str or sympy expression) for a single expression 
                          or an instance of the :py:class:`genepy.Gaussian` class
    
    Methods:
    --------
  """
  def __init__(self,hs_expression) -> None:
    self.data = hs_expression
    return
  
  def __str__(self) -> str:
    return self.data.__str__()
  
  def sprint_option_gaussian(self, model_name: str):
    """
    sprint_option_gaussian(self,model_name:str)
    Return a string formatted for `pTatin3d`_ input file using `PETSc`_ options format.

    :param str model_name: name of the model

    :return: string with the options
    """
    prefix = "heat_source"
    s  = f"########### Initial heat source ###########\n"
    s += self.data.sprint_option(model_name, prefix)
    return s
  
  def sprint_option(self, model_name: str):
    """
    sprint_option(self,model_name:str)
    Return a string formatted for `pTatin3d`_ input file using `PETSc`_ options format.

    :param str model_name: name of the model

    :return: string with the options
    """
    if isinstance(self.data, Gaussian):
      return self.sprint_option_gaussian(model_name)
    else:
      prefix = "heat_source"
      s  = f"########### Initial heat source ###########\n"
      s += f"-{model_name}_{prefix}_nhs 1\n"
      s += f"-{model_name}_{prefix}_expression {str(self.data).replace(' ','')}\n"
      return s