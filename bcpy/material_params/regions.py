from bcpy import MaterialConstants

# This file contains the Region class which is used to store the material parameters for a given region of the model.
class Region:
  def __init__(self, region: int, 
         density: MaterialConstants, 
         plasticity: MaterialConstants,
         softening: MaterialConstants, 
         viscosity: MaterialConstants,
         energy: MaterialConstants=None) -> None:
    
    self.region = region
    self.material_parameters = [density, plasticity, softening, viscosity]
    if energy is not None: self.material_parameters.append(energy)
    # enforce region number for all material parameters
    for m in self.material_parameters:
      m.region = region

  def sprint_option(self) -> str:
    s = f"###### Region {self.region} ######\n"
    for m in self.material_parameters:
      s += m.sprint_option()
    return s
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}: {self.region}\n'
    for m in self.material_parameters:
      s += str(m)
    return s
  
class ModelRegions:
  def __init__(self, model_name:str, regions:list) -> None:
    self.model_name = model_name
    self.regions    = regions

    for r in self.regions:
      for m in r.material_parameters:
        m.model_name = self.model_name
  
  def sprint_option(self) -> str:
    prefix = "regions"
    nregions = len(self.regions)
    s  = "########### Material parameters ###########\n"
    s += f"-{self.model_name}_{prefix}_nregions {nregions}\n"
    s += f"-{self.model_name}_{prefix}_list "
    for r in range(nregions-1):
      s += f"{self.regions[r].region},"
    s += f"{self.regions[nregions-1].region}\n"

    for r in self.regions:
      s += r.sprint_option()
    return s
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    for r in self.regions:
      s += str(r)
    return s
  