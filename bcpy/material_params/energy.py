from bcpy import MaterialConstants

class Energy(MaterialConstants):
  def __init__(self, heat_source:float, conductivity:float, Cp:float=800.0, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heatcapacity = Cp
    self.heat_source  = heat_source
    self.conductivity = conductivity
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:        {self.region}\n'
    s += f'\tHeat capacity: {self.heatcapacity}\n'
    s += f'\tHeat source:   {self.heat_source}\n'
    s += f'\tConductivity:  {self.conductivity}\n'
    return s