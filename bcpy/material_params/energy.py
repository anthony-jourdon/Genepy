#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  bc-pre-processing
#  filename: energy.py
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

class Energy(MaterialConstants):
  """
  .. py:class:: Energy(heat_source:float, conductivity:float, Cp:float=800.0, model_name:str="model_GENE3D", region:int=0)

    Class to define the energy parameters of a region of the model.

    :param float heat_source: Heat source value in :math:`W.m^{-3}`.
    :param float conductivity: Conductivity value in :math:`W.m^{-1}.K^{-1}`.
    :param float Cp: Heat capacity value in :math:`J.kg^{-1}.K^{-1}`. Default is :math:`800` :math:`J.kg^{-1}.K^{-1}`.
    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.

    Attributes
    ----------

    .. py:attribute:: heatcapacity
      :type: float

        Heat capacity value in :math:`J.kg^{-1}.K^{-1}`. Default is :math:`800` :math:`J.kg^{-1}.K^{-1}`.

    .. py:attribute:: heat_source
      :type: float

        Heat source value in :math:`W.m^{-3}`.

    .. py:attribute:: conductivity
      :type: float

        Conductivity value in :math:`W.m^{-1}.K^{-1}`.
  """
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
