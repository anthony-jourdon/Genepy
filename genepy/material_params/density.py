#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: density.py
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

from genepy import MaterialConstants

class Density(MaterialConstants):
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self,model_name,region)

class DensityConstant(Density):
  """
  .. py:class:: DensityConstant(density:float, model_name:str="model_GENE3D", region:int=0)

    Class to apply a constant density to a region of the model.

    :param float density: Constant density value in kg.m\ :sup:`-3` to be applied to the region.
    :param str model_name: Name of the model to which the density is applied. Default is ``"model_GENE3D"``.
    :param int region: Region number to which the density is applied. Default is 0.

    .. note::
      The density is given by:
        :math:`\\rho(p,T) = \\rho_0`

    Attributes
    ----------
    
    .. py:attribute:: density
      :type: float

        Constant density value (:math:`\\rho_0`) in kg.m\ :sup:`-3` to be applied to the region.
    
    .. py:attribute:: model_name
      :type: str

        Name of the model to which the density is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the density is applied. Default is 0.
  """
  def __init__(self, density:float, model_name:str="model_GENE3D", region:int=0) -> None:
    self.density_type = 0
    self.density      = density
    Density.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:  {self.region}\n'
    s += f'\tDensity: {self.density}\n'
    return s
  
class DensityBoussinesq(Density):
  """
  .. py:class:: DensityBoussinesq(density:float, thermal_expansion:float=0.0, compressibility:float=0.0, model_name:str="model_GENE3D", region:int=0)

    Class to apply a density using the Boussinesq approximation to a region of the model.

    :param float density: Reference density value in kg.m\ :sup:`-3` (:math:`\\rho_0`).
    :param float thermal_expansion: Thermal expansion coefficient in K\ :sup:`-1` (:math:`\\alpha`).
    :param float compressibility: Compressibility coefficient in Pa\ :sup:`-1` (:math:`\\beta`).
    :param str model_name: Name of the model to which the density is applied. Default is "model_GENE3D".
    :param int region: Region number to which the density is applied. Default is 0.

    .. note::
      The density is given by:
        :math:`\\rho(p,T) = \\rho_0 (1 - \\alpha T + \\beta p)`

    Attributes
    ----------

    .. py:attribute:: density
      :type: float

        Reference density value in kg.m\ :sup:`-3` (:math:`\\rho_0`).

    .. py:attribute:: thermal_expansion
      :type: float

        Thermal expansion coefficient in K\ :sup:`-1` (:math:`\\alpha`).

    .. py:attribute:: compressibility
      :type: float

        Compressibility coefficient in Pa\ :sup:`-1` (:math:`\\beta`).

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the density is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the density is applied. Default is 0.
  """
  def __init__(self, density:float, thermal_expansion:float=0.0, compressibility:float=0.0, model_name:str="model_GENE3D", region:int=0) -> None:
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
  """
  .. py:class:: DensityTable(density:float, map:str, model_name:str="model_GENE3D", region:int=0)

    Class to apply a density read from a thermodynamic table to a region of the model.

    :param float density: Density value in kg.m\ :sup:`-3` in case the map does not contain the :math:`p,T` conditions.
    :param str map: Name of the file containing the density map.
    :param str model_name: Name of the model to which the density is applied. Default is "model_GENE3D".
    :param int region: Region number to which the density is applied. Default is 0.

    Attributes
    ----------

    .. py:attribute:: density
      :type: float

        Density value in kg.m\ :sup:`-3` in case the map does not contain the :math:`p,T` conditions.

    .. py:attribute:: map
      :type: str

        Name of the file containing the density map.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the density is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the density is applied. Default is 0.
  """
  def __init__(self, density:float, map:str, model_name:str="model_GENE3D", region:int=0) -> None:
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
