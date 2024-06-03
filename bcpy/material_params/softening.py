#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: softening.py
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

from bcpy import MaterialConstants

class Softening(MaterialConstants):
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self,model_name,region)

class SofteningNone(Softening):
  """
  .. py:class:: SofteningNone(model_name:str="model_GENE3D", region:int=0)

    Class to apply no softening to a region of the model.
    The given region will not have any softening.

    :param str model_name: Name of the model to which the softening is applied. Default is "model_GENE3D".
    :param int region: Region number to which the softening is applied. Default is 0.

  """
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    self.softening_type = 0
    Softening.__init__(self,model_name,region)

class SofteningLinear(Softening):
  """
  .. py:class:: SofteningLinear(strain_min:float, strain_max:float, model_name:str='model_GENE3D', region:int=0)

    Class to apply linear softening to a region of the model.
    Only applies to the plastic deformation.

    :param float strain_min: Minimum strain value (:math:`\\epsilon_{min}`).
    :param float strain_max: Maximum strain value (:math:`\\epsilon_{max}`).
    :param str model_name: Name of the model to which the softening is applied. Default is "model_GENE3D".
    :param int region: Region number to which the softening is applied. Default is 0.
  
    .. note:: 
      The linear softening is defined by the minimum and maximum plastic strain values and influences the friction and/or the cohesion.
      For :math:`t` the variable to be softened, the softening is given by:
      
      .. math::
        t = t_0 - (t_0 - t_{\\infty}) \\frac{\\epsilon_p - \\epsilon_{min}}{\\epsilon_{max} - \\epsilon_{min}}
      
      where :math:`t_0` is the initial value of the variable, :math:`t_{\\infty}` is the final value of the variable, and :math:`\\epsilon_p` is the plastic strain.
  
    Attributes
    ----------

    .. py:attribute:: eps_min
      :type: float
        
        Minimum strain value (:math:`\\epsilon_{min}`).

    .. py:attribute:: eps_max
      :type: float
        
        Maximum strain value (:math:`\\epsilon_{max}`).

    .. py:attribute:: model_name
      :type: str
        
        Name of the model to which the softening is applied.

    .. py:attribute:: region
      :type: int
        
        Region number to which the softening is applied.
  
  """
  def __init__(self, strain_min:float, strain_max:float, model_name:str="model_GENE3D", region:int=0) -> None:
    self.softening_type = 1
    self.eps_min        = strain_min
    self.eps_max        = strain_max
    Softening.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'\tRegion:    {self.region}\n'
    s += f'\tStrain min: {self.eps_min}\n'
    s += f'\tStrain max: {self.eps_max}\n'
    return s

class SofteningExponential(Softening):
  """
  .. py:class:: SofteningExponential(strain_min:float, strain_fold:float, model_name:str='model_GENE3D', region:int=0)

    Class to apply exponential softening to a region of the model.
    Only applies to the plastic deformation.

    :param float strain_min: Minimum strain value (:math:`\\epsilon_{min}`).
    :param float strain_fold: Fold strain value (:math:`\\epsilon_{fold}`).
    :param str model_name: Name of the model to which the softening is applied. Default is "model_GENE3D".
    :param int region: Region number to which the softening is applied. Default is 0.

    
    .. note:: 
      The exponential softening is defined by the minimum and a fold plastic strain values and influences the friction and/or the cohesion.
      For :math:`t` the variable to be softened, the softening is given by:
        
      .. math:: 
        t = t_0 - (t_0 - t_{\\infty}) \\exp \\left(-\\frac{\\epsilon_p}{\\epsilon_{fold}} \\right)
      
      where :math:`t_0` is the initial value of the variable, :math:`t_{\\infty}` is the final value of the variable, and :math:`\\epsilon_p` is the plastic strain.

    Attributes
    ----------

    .. py:attribute:: eps_min
      :type: float
        
        Minimum strain value (\\epsilon_{min}).

    .. py:attribute:: eps_fold
      :type: float
        
        Fold strain value (\\epsilon_{fold}).

    .. py:attribute:: model_name
      :type: str
        
        Name of the model to which the softening is applied.

    .. py:attribute:: region
      :type: int

        Region number to which the softening is applied.
  """
  def __init__(self, strain_min:float, strain_fold:float, model_name:str="model_GENE3D", region:int=0) -> None:
    self.softening_type = 2
    self.eps_min        = strain_min
    self.eps_fold       = strain_fold
    Softening.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'\tRegion:      {self.region}\n'
    s += f'\tStrain min:  {self.eps_min}\n'
    s += f'\tStrain fold: {self.eps_fold}\n'
    return s
