from bcpy import MaterialConstants

class Softening(MaterialConstants):
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self,model_name,region)

class SofteningNone(Softening):
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    self.softening_type = 0
    Softening.__init__(self,model_name,region)

class SofteningLinear(Softening):
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
