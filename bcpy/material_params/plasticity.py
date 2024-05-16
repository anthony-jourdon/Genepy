from bcpy import MaterialConstants

class Plasticity(MaterialConstants):
  def __init__(self, model_name:str, region:int) -> None:
    MaterialConstants.__init__(self,model_name,region)

class PlasticNone(Plasticity):
  def __init__(self, model_name: str, region: int) -> None:
    self.plastic_type = 0
    Plasticity.__init__(self,model_name,region)

class PlasticMises(Plasticity):
  def __init__(self, model_name:str, region:int, yield_stress:float, yield_stress_inf:float) -> None:
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
  def __init__(self, model_name:str, region:int, friction:float, friction_inf:float, cohesion:float, cohesion_inf:float, tension_cutoff:float, highstress_cutoff:float) -> None:
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