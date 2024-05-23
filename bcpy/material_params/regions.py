from bcpy import MaterialConstants
from bcpy import PlasticNone
from bcpy import DensityConstant
from bcpy import ViscosityConstant
from bcpy import SofteningNone

# This file contains the Region class which is used to store the material parameters for a given region of the model.
class Region:
  def __init__(self, region: int, 
         density: MaterialConstants=None, 
         viscosity: MaterialConstants=None,
         plasticity: MaterialConstants=None,
         softening: MaterialConstants=None, 
         energy: MaterialConstants=None) -> None:
    
    self.region = region

    # if user does not provide anything, set default values
    if density is None:    density = DensityConstant(3300.0)
    if viscosity is None:  viscosity = ViscosityConstant(1.0e22)
    if plasticity is None: plasticity = PlasticNone()
    if softening is None:  softening = SofteningNone()

    self.material_parameters = [density, viscosity, plasticity, softening]
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
  def __init__(self, regions:list[MaterialConstants], model_name:str="model_GENE3D") -> None:
    self.model_name = model_name
    self.regions    = regions

    for r in self.regions:
      for m in r.material_parameters:
        m.model_name = self.model_name
  
  def add_region(self, region:Region) -> None:
    for m in region.material_parameters:
      m.model_name = self.model_name
    self.regions.append(region)
    return
  
  def sprint_option(self) -> str:
    prefix = "regions"
    nregions = len(self.regions)
    s  = "########### Material parameters ###########\n"
    s += f"-{self.model_name}_mesh_file path_to_file\n"
    s += f"-{self.model_name}_{prefix}_file path_to_file\n"
    s += f"-{self.model_name}_{prefix}_nregions {nregions}\n"
    s += f"-{self.model_name}_{prefix}_list "
    for r in range(nregions-1):
      s += f"{self.regions[r].region},"
    s += f"{self.regions[nregions-1].region}\n"
    s += "# Method to locate material points in gmsh mesh\n"
    s += "# Brute force: 0, Partitioned box: 1\n"
    s += f"-{self.model_name}_mesh_point_location_method 1\n"
    for r in self.regions:
      s += r.sprint_option()
    return s
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    for r in self.regions:
      s += str(r)
    return s
  