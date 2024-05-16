from bcpy import MaterialConstants

class Density(MaterialConstants):
  def __init__(self, model_name:str, region:int) -> None:
    MaterialConstants.__init__(self,model_name,region)

class DensityConstant(Density):
  def __init__(self, model_name:str, region:int, density:float) -> None:
    self.density_type = 0
    self.density      = density
    Density.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:  {self.region}\n'
    s += f'\tDensity: {self.density}\n'
    return s
  
class DensityBoussinesq(Density):
  def __init__(self, model_name: str, region: int, density:float, thermal_expansion:float=0.0, compressibility:float=0.0) -> None:
    self.density_type     = 1
    self.density          = density
    self.thermalexpension = thermal_expansion
    self.compressibility  = compressibility
    Density.__init__(self,model_name, region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'rho(p,T) = rho_0 * (1 - alpha * T + beta * p)\n'
    s += f'\tRegion:            {self.region}\n'
    s += f'\tDensity:           {self.density}\n'
    s += f'\tThermal expansion: {self.thermalexpension}\n'
    s += f'\tCompressibility:   {self.compressibility}\n'
    return s

class DensityTable(Density):
  def __init__(self, model_name: str, region: int, density:float, map:str) -> None:
    self.density_type = 2
    self.density      = density
    self.map          = map
    Density.__init__(self,model_name, region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'Density read from thermodynamic table\n'
    s += f'\tRegion:  {self.region}\n'
    s += f'\tDensity: {self.density}\n'
    s += f'\tMap:     {self.map}\n'
    return s

def test():
  region = 0
  dt = DensityConstant("model_GENE3D",region,2700.0)
  region = 1
  db = DensityBoussinesq("model_GENE3D",region,2900.0,3e-5,1e-11)
  region = 2
  dtb = DensityTable("model_GENE3D",region,3300.0,"density_map.txt")

  opt  = dt.sprint_option()
  opt += db.sprint_option()
  opt += dtb.sprint_option()
  print(opt)

if __name__ == "__main__":
  test()