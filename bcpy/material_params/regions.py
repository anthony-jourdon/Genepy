#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: regions.py
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
from bcpy import PlasticNone
from bcpy import DensityConstant
from bcpy import ViscosityConstant
from bcpy import SofteningNone

# This file contains the Region class which is used to store the material parameters for a given region of the model.
class Region:
  """
  .. py:class:: Region(region:int, density:MaterialConstants=None, viscosity:MaterialConstants=None, plasticity:MaterialConstants=None, softening:MaterialConstants=None, energy:MaterialConstants=None)

    Class to define the material parameters of a region of the model.
    The region number is enforced for all material parameters.

    :param int region: Region number to which the material parameters are applied.
    :param MaterialConstants density: instance of :ref:`Density <density>` class of the region.
    :param MaterialConstants viscosity: instance of :ref:`Viscosity <viscosity>` class of the region.
    :param MaterialConstants plasticity: instance of :ref:`Plasticity <plasticity>` class of the region.
    :param MaterialConstants softening: instance of :ref:`Softening <softening>` class of the region.
    :param MaterialConstants energy: instance of :ref:`Energy <energy>` class of the region.

    .. note::
      If the user does not provide any material parameter, the default values are set to:
      
      - :py:class:`Constant Density <bcpy.material_params.density.DensityConstant>` : 3300 kg/m\ :sup:`3`
      - :py:class:`Constant Viscosity <bcpy.material_params.viscosity.ViscosityConstant>` : 10\ :sup:`22` Pa.s
      - :py:class:`No Plasticity <bcpy.material_params.plasticity.PlasticNone>` 
      - :py:class:`No Softening <bcpy.material_params.softening.SofteningNone>` 
    
    Example
    -------

    Default values for all material parameters
    ..........................................

    >>> import bcpy as bp
    >>> region = bp.Region(1)

    Custom values for material parameters
    .....................................

    .. code-block:: python

      import bcpy as bp

      region1 = bp.Region(1,                                           # region tag
                          bp.DensityBoussinesq(2700.0,3.0e-5,1.0e-11), # density
                          bp.ViscosityArrhenius2("Quartzite"),         # viscosity  (values from the database using rock name)
                          bp.SofteningLinear(0.0,0.5),                 # softening
                          bp.PlasticDruckerPrager(),                   # plasticity (default values, can be modified using the corresponding parameters)
                          bp.Energy(1.0e-6,2.7)),                      # energy

  """
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
    """
    sprint_option(self)
    :return: String containing all the options for all material parameters of the given region.
    :rtype: str
    """
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
  """
  .. py:class:: ModelRegions(regions:list[Region], mesh_file:str="path_to_file", region_file:str="path_to_file", model_name:str="model_GENE3D")

    Class to define the material parameters of a model.
    The user can provide a list of instances of class :py:class:`Region` with their material parameters.
    The model name is enforced for all material parameters.
    
    :param list[Region] regions: List of instances of class :py:class:`Region` with their material parameters.
    :param str mesh_file: (**Optional**) path to the binary mesh file.
    :param str region_file: (**Optional**) Path to the binary region file.
    :param str model_name: Name of the model. Default is "model_GENE3D".
  """
  def __init__(self, regions:list[MaterialConstants], mesh_file:str="path_to_file", region_file:str="path_to_file", model_name:str="model_GENE3D") -> None:
    self.model_name  = model_name
    self.mesh_file   = mesh_file
    self.region_file = region_file
    self.regions     = regions

    for r in self.regions:
      for m in r.material_parameters:
        m.model_name = self.model_name
  
  def add_region(self, region:Region) -> None:
    """
    add_region(self, region:Region)
    Adds a region to the list of regions after the instance of the class has been created.
    
    :param Region region: Instance of class py:class:`Region` with its material parameters.
    """
    for m in region.material_parameters:
      m.model_name = self.model_name
    self.regions.append(region)
    return
  
  def sprint_option(self) -> str:
    """
    sprint_option(self)
    :return: String containing all the options for all regions of the model.
    :rtype: str
    """
    prefix = "regions"
    nregions = len(self.regions)
    s  = "########### Material parameters ###########\n"
    s += f"-{self.model_name}_mesh_file {self.mesh_file}\n"
    s += f"-{self.model_name}_{prefix}_file {self.region_file}\n"
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
  
