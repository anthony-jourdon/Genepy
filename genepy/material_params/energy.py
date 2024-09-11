#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: energy.py
#
#  This file is part of Genepy.
#
#  Genepy is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  Genepy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with Genepy. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

from genepy import MaterialConstants

"""
  ENERGYSOURCE_NONE = 0,
  ENERGYSOURCE_USE_MATERIALPOINT_VALUE,
  ENERGYSOURCE_CONSTANT ,
  ENERGYSOURCE_SHEAR_HEATING,
  ENERGYSOURCE_DECAY,
  ENERGYSOURCE_ADIABATIC,
  ENERGYSOURCE_ADIABATIC_ADVECTION
"""

class EnergySourceNone(MaterialConstants):
  """
  .. py:class:: EnergySourceNone(model_name:str="model_GENE3D", region:int=0)

    Class to define the energy source type of a region of the model as none.
    Default behavior of `pTatin3d`_.

    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.
  """
  def __init__(self, model_name: str = "model_GENE3D", region: int = 0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heat_source_type = 0

  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'\tRegion: {self.region}\n'
    return s

class EnergySourceMaterialPointValue(MaterialConstants):
  """
  .. py:class:: EnergySourceMaterialPointValue(model_name:str="model_GENE3D", region:int=0)

    Class to define the energy source type of a region of the model as a defined on material points.
    This class only sets the type of heat source to be used. 
    The actual values of the heat source to apply on the material points should be defined in the
    :py:class:`genepy.InitialHeatSource` class with an expression or an instance of the  
    :py:class:`genepy.Gaussian` class describing the heat source spatial repartition.

    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.
  """
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heat_source_type = 1
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'\tRegion: {self.region}\n'
    return s

class EnergySourceConstant(MaterialConstants):
  """
  .. py:class:: EnergySourceConstant(heat_source:float, model_name:str="model_GENE3D", region:int=0)

    Class to define the energy source type of a region of the model as a constant value.
    The heat source value should be given in W.m\\ :sup:`-3`.

    :param float heat_source: Heat source value in W.m\\ :sup:`-3`.
    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.
  """
  def __init__(self, heat_source:float, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heat_source_type = 2
    self.heat_source      = heat_source
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'\tRegion:      {self.region}\n'
    s += f'\tHeat source: {self.heat_source}\n'
    return s

class EnergySourceShearHeating(MaterialConstants):
  """
  .. py:class:: EnergySourceShearHeating(model_name:str="model_GENE3D", region:int=0)

    Class to define the energy source type of a region of the model as shear heating.
    The source derived from shear heating is defined as:

    .. math:: 
      H &= \\boldsymbol{\\sigma} : \\boldsymbol{\\varepsilon} \\\\
        &= 2 \\eta \\boldsymbol{\\varepsilon} : \\boldsymbol{\\varepsilon} \\\\
        &= 2 \\eta \\sum_{ij} \\varepsilon_{ij} \\varepsilon_{ij}

    with :math:`\\eta` the viscosity and :math:`\\boldsymbol{\\varepsilon}` the strain rate tensor.

    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.
  """
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heat_source_type = 3

  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion: {self.region}\n'
    return s

class EnergySourceDecay(MaterialConstants):
  """
  .. py:class:: EnergySourceDecay(heat_source:float, half_life:float, model_name:str="model_GENE3D", region:int=0)

    Class to define the energy source type of a region of the model as a decaying heat source such that

    .. math:: 
      H_0 e^{-\\lambda t}
    
    where :math:`H_0` is the reference heat source value in W.m\\ :sup:`-3`, 
    :math:`\\lambda` is the decay constant in s\\ :sup:`-1` and :math:`t` is the time.

    :param float heat_source: Reference heat source value :math:`H_0` in W.m\\ :sup:`-3`.
    :param float half_life: Decay constant in s\\ :sup:`-1`
    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.
  """
  def __init__(self, heat_source:float, decay:float, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heat_source_type      = 4
    self.heat_source_ref       = heat_source
    self.heat_source_half_life = decay
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:          {self.region}\n'
    s += f'\tHeat source ref: {self.heat_source_ref}\n'
    s += f'\tHalf life:       {self.heat_source_half_life}\n'
    return s

class EnergySource(MaterialConstants):
  """
  .. py:class:: EnergySource(*args, model_name:str="model_GENE3D", region:int=0)

    Class to define the energy source type of a region of the model as a combination of different heat sources.
    Those sources will be summed in the model such that the total heat source at a point in space and time is defined by

    .. math::
      H(\\mathbf x, t) = \\sum_i H_i(\\mathbf x, t)

    where :math:`H_i(\\mathbf x, t)` is the heat source value of the :math:`i`-th source type. 
    This way of defining the heat source allows to combine different types of heat sources in the model.

    Example:
    --------

    To define a region with a combination of heat sources defined on 
    :py:class:`material points <genepy.EnergySourceMaterialPointValue>`, 
    plus a :py:class:`constant <genepy.EnergySourceConstant>` of 1.5 :math:`\\mu`\\ W.m\\ :sup:`-3` 
    and a :py:class:`shear heating <genepy.EnergySourceShearHeating>` source, one can use:

    .. code-block:: python

      import genepy as gp

      gp.EnergySource(gp.EnergySourceMaterialPointValue(),
                      gp.EnergySourceConstant(1.5e-6),
                      gp.EnergySourceShearHeating())

    :param args: Heat sources to combine for the given region, can pass multiple heat sources. 
                Can only be of type: 
                 
                - :py:class:`genepy.EnergySourceMaterialPointValue`
                - :py:class:`genepy.EnergySourceConstant`
                - :py:class:`genepy.EnergySourceShearHeating`
                - :py:class:`genepy.EnergySourceDecay`

    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.
  """
  def __init__(self, *args, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heat_sources = args
    for arg in self.heat_sources:
      if not isinstance(arg, (EnergySourceNone, EnergySourceMaterialPointValue, EnergySourceConstant, EnergySourceShearHeating, EnergySourceDecay)):
        raise ValueError(f"arg, expected one of {EnergySourceMaterialPointValue}, {EnergySourceConstant}, {EnergySourceShearHeating}, {EnergySourceDecay}, got {type(arg)}")
    return

  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    for hs in self.heat_sources:
      s += hs.__str__()
    return s

  def sprint_option(self):
    s = f"-{self.model_name}_heat_source_type_{self.region} "
    for i in range(len(self.heat_sources)-1):
      s += f"{self.heat_sources[i].heat_source_type},"
    s += f"{self.heat_sources[-1].heat_source_type}\n"

    for hs in self.heat_sources:
      s += hs.sprint_option()
    return s

class Energy(MaterialConstants):
  """
  .. py:class:: Energy(heat_source:EnergySource, conductivity:float, Cp:float=800.0, model_name:str="model_GENE3D", region:int=0)

    Class to define the energy parameters of a region of the model.

    :param EnergySource heat_source: Instance of the :py:class:`genepy.EnergySource` class describing the heat source types.
    :param float conductivity: Conductivity value in W.m\\ :sup:`-1`.K\\ :sup:`-1`.
    :param float Cp: Heat capacity value in J.kg\\ :sup:`-1`.K\\ :sup:`-1`. Default is :math:`800` J.kg\\ :sup:`-1`.K\\ :sup:`-1`.
    :param str model_name: Name of the model to which the energy parameters are applied. Default is "model_GENE3D".
    :param int region: Region number to which the energy parameters are applied. Default is 0.

    Attributes
    ----------

    .. py:attribute:: heatcapacity
      :type: float

        Heat capacity value in J.kg\\ :sup:`-1`.K\\ :sup:`-1`. Default is :math:`800` J.kg\\ :sup:`-1`.K\\ :sup:`-1`.

    .. py:attribute:: heat_source
      :type: float

        Instance of the :py:class:`genepy.EnergySource` class describing the heat source types.

    .. py:attribute:: conductivity
      :type: float

        Conductivity value in W.m\\ :sup:`-1`.K\\ :sup:`-1`.
  """
  def __init__(self, heat_source:EnergySource, conductivity:float, Cp:float=800.0, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self, model_name, region)
    self.heatcapacity = Cp
    self.heat_source  = heat_source
    self.conductivity = conductivity
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:        {self.region}\n'
    s += f'\tHeat capacity: {self.heatcapacity}\n'
    s += f'\tConductivity:  {self.conductivity}\n'
    s += self.heat_source.__str__()
    return s
  
  def sprint_option(self):
    self.heat_source.region = self.region
    for hs in self.heat_source.heat_sources:
      hs.region = self.region

    s = f"-{self.model_name}_heatcapacity_{self.region} {self.heatcapacity}\n"
    s += f"-{self.model_name}_conductivity_{self.region} {self.conductivity}\n"
    s += self.heat_source.sprint_option()
    return s