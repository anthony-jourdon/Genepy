import numpy as np
from bcpy import MaterialConstants

class Plasticity(MaterialConstants):
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self,model_name,region)

class PlasticNone(Plasticity):
  """
  .. py:class:: PlasticNone(model_name:str="model_GENE3D", region:int=0)

    Class to apply no plasticity to a region of the model.
    The given region will not have any plastic behavior.

    :param str model_name: Name of the model to which the plasticity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the plasticity is applied. Default is 0.

    Attributes
    ----------

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the plasticity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the plasticity is applied. Default is 0
  """
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    self.plastic_type = 0
    Plasticity.__init__(self,model_name,region)

class PlasticMises(Plasticity):
  """
  .. py:class:: PlasticMises(yield_stress:float=5.e7, yield_stress_inf:float=5.e6, model_name:str="model_GENE3D", region:int=0)

    Class to apply von Mises plasticity to a region of the model.

    :param float yield_stress: Yield stress value in Pa (C). Default is :math:`50` MPa.
    :param float yield_stress_inf: Yield stress value in Pa after softening (if any). Default is :math:`5` MPa.
    :param str model_name: Name of the model to which the plasticity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the plasticity is applied. Default is 0.

    .. note:: 
      The yield stress is given by :math:`\\sigma_y = C`

    Attributes
    ----------

    .. py:attribute:: yield_stress
      :type: float

        Yield stress value in Pa (C). Default is :math:`50` MPa.

    .. py:attribute:: yield_stress_inf
      :type: float

        Yield stress value in Pa after softening (if any). Default is :math:`5` MPa.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the plasticity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the plasticity is applied. Default is 0
  """
  def __init__(self, yield_stress:float=5.e7, yield_stress_inf:float=5.e6,
               model_name:str="model_GENE3D", region:int=0) -> None:
    self.plastic_type     = 1
    self.yield_stress     = yield_stress
    self.yield_stress_inf = yield_stress_inf
    Plasticity.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'\tRegion:          {self.region}\n'
    s += f'\tYield stress:    {self.yield_stress}\n'
    s += f'\tYield stress inf:{self.yield_stress_inf}\n'
    return s

class PlasticDruckerPrager(Plasticity):
  """
  .. py:class:: PlasticDruckerPrager(friction:float=np.deg2rad(30.), friction_inf:float=np.deg2rad(5.), cohesion:float=2.e7, cohesion_inf:float=5.e6, tension_cutoff:float=1.e7, highstress_cutoff:float=4.e8, model_name:str="model_GENE3D", region:int=0)

    Class to apply Drucker-Prager plasticity to a region of the model.

    .. note:: 
      The yield stress is given by 
        :math:`\\sigma_y = C \\cos(\\phi) + p \\sin(\\phi)`

      with :math:`\\phi` the friction angle in radian and :math:`C` the cohesion in Pa.

    :param float friction: Friction angle value in radian (:math:`\\phi`). Default is :math:`30` degrees.
    :param float friction_inf: Friction angle value in radian after softening (if any). Default is :math:`5` degrees.
    :param float cohesion: Cohesion value in Pa (:math:`C`). Default is :math:`20` MPa.
    :param float cohesion_inf: Cohesion value in Pa after softening (if any). Default is :math:`5` MPa.
    :param float tension_cutoff: Tension cutoff value in Pa (if :math:`\\sigma_y < \\text{tension_cutoff} : \\sigma_y = \\text{tension_cutoff}`). Default is :math:`10` MPa.
    :param float highstress_cutoff: High stress cutoff value in Pa (if :math:`\\sigma_y > \\text{highstress_cutoff} : \\sigma_y = \\text{highstress_cutoff}`). Default is :math:`400` MPa.
    :param str model_name: Name of the model to which the plasticity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the plasticity is applied. Default is 0.

    Attributes
    ----------

    .. py:attribute:: friction
      :type: float

        Friction angle value in radian (:math:`\\phi`). Default is :math:`30` degrees.

    .. py:attribute:: friction_inf
      :type: float

        Friction angle value in radian after softening (if any). Default is :math:`5` degrees.

    .. py:attribute:: cohesion
      :type: float

        Cohesion value in Pa (:math:`C`). Default is :math:`20` MPa.

    .. py:attribute:: cohesion_inf
      :type: float

        Cohesion value in Pa after softening (if any). Default is :math:`5` MPa.

    .. py:attribute:: tension_cutoff
      :type: float

        Tension cutoff value in Pa (if :math:`\\sigma_y < \\text{tension_cutoff} : \\sigma_y = \\text{tension_cutoff}`). Default is :math:`10` MPa.

    .. py:attribute:: highstress_cutoff
      :type: float

        High stress cutoff value in Pa (if :math:`\\sigma_y > \\text{highstress_cutoff} : \\sigma_y = \\text{highstress_cutoff}`). Default is :math:`400` MPa.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the plasticity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the plasticity is applied. Default is 0
  """
  def __init__(self, friction:float=np.deg2rad(30.), friction_inf:float=np.deg2rad(5.),
               cohesion:float=2.e7, cohesion_inf:float=5.e6, 
               tension_cutoff:float=1.e7, highstress_cutoff:float=4.e8,
               model_name:str="model_GENE3D", region:int=0) -> None:
    self.plastic_type = 2
    self.friction     = friction
    self.friction_inf = friction_inf
    self.cohesion     = cohesion
    self.cohesion_inf = cohesion_inf
    self.tens_cutoff  = tension_cutoff
    self.hst_cutoff   = highstress_cutoff
    Plasticity.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'sigma_y = C * cos(phi) + p * sin(phi)\n'
    s += f'\tRegion:             {self.region}\n'
    s += f'\tFriction:           {self.friction}\n'
    s += f'\tFriction inf:       {self.friction_inf}\n'
    s += f'\tCohesion:           {self.cohesion}\n'
    s += f'\tCohesion inf:       {self.cohesion_inf}\n'
    s += f'\tTension cutoff:     {self.tens_cutoff}\n'
    s += f'\tHigh-stress cutoff: {self.hst_cutoff}\n'
    return s