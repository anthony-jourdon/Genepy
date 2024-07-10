from genepy.initial_conditions.gaussian import Gaussian

class InitialPlasticStrain:
  """
  .. py:class:: InitialPlasticStrain(strain_expressions)

    Class to generate options for the initial plastic strain. 

    :param strain_expressions: instance of the :py:class:`genepy.Gaussian` class describing the initial plastic strain
    
    Methods:
    --------
    """
  def __init__(self,strain_expressions:Gaussian) -> None:
    self.data = strain_expressions
    return
  
  def __str__(self) -> str:
    return self.data.__str__()

  def sprint_option(self, model_name: str):
    """
    sprint_option(self,model_name:str)
    Return a string formatted for `pTatin3d`_ input file using `PETSc`_ options format.

    :param str model_name: name of the model

    :return: string with the options
    """
    prefix = "wz"
    s  = f"########### Initial plastic strain for weak zone ###########\n"
    s += self.data.sprint_option(model_name, prefix)
    return s