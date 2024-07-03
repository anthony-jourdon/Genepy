#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  Genepy
#  filename: viscosity.py
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

import os
import json
from genepy import MaterialConstants

class Viscosity(MaterialConstants):
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    MaterialConstants.__init__(self,model_name,region)

class ViscosityConstant(Viscosity):
  """
  .. py:class:: ViscosityConstant(viscosity:float, model_name:str="model_GENE3D", region:int=0)

    Class to apply a constant viscosity to a region of the model.

    :param float viscosity: Constant viscosity value in Pa.s to be applied to the region.
    :param str model_name: Name of the model to which the viscosity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the viscosity is applied. Default is 0.

    .. note::
      The viscosity is given by:
        :math:`\\eta = \\eta_0`

    Attributes
    ----------

    .. py:attribute:: viscosity
      :type: float

        Constant viscosity value (:math:`\\eta_0`) in Pa.s to be applied to the region.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the viscosity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the viscosity is applied. Default is 0.
  """
  def __init__(self, viscosity:float, model_name:str="model_GENE3D", region:int=0) -> None:
    self.visc_type = 0
    self.eta0           = viscosity
    Viscosity.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:    {self.region}\n'
    s += f'\tViscosity: {self.eta0}\n'
    return s
  
class ViscosityFrankK(Viscosity):
  """
  .. py:class:: ViscosityFrankK(eta0:float, exponent:float, model_name:str="model_GENE3D", region:int=0)

    Class to apply a Frank-Kamenetskii viscosity to a region of the model.

    :param float eta0: Viscosity value at the reference temperature (:math:`\\eta_0`) in Pa.s.
    :param float exponent: Exponent of the viscosity law (:math:`\\theta`).
    :param str model_name: Name of the model to which the viscosity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the viscosity is applied. Default is 0.

    .. note::
      The viscosity is given by:
        :math:`\\eta(T) = \\eta_0 \\exp(-\\theta T)`

    Attributes
    ----------

    .. py:attribute:: eta0
      :type: float

        Viscosity value at the reference temperature (:math:`\\eta_0`) in Pa.s.

    .. py:attribute:: exponent
      :type: float

        Exponent of the viscosity law (:math:`\\theta`).

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the viscosity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the viscosity is applied. Default is 0.
  """
  def __init__(self, eta0:float, exponent:float, model_name:str="model_GENE3D", region:int=0) -> None:
    self.visc_type = 1
    self.eta0           = eta0
    self.theta          = exponent
    Viscosity.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'eta(T) = eta_0 * exp(-theta * T)\n'
    s += f'\tRegion:   {self.region}\n'
    s += f'\tViscosity:{self.eta0}\n'
    s += f'\tExponent: {self.theta}\n'
    return s
  
class ViscosityZ(Viscosity):
  """
  .. py:class:: ViscosityZ(eta0:float, zeta:float, zref:float, model_name:str="model_GENE3D", region:int=0)

    Class to apply a depth dependant viscosity to a region of the model.

    :param float eta0: Viscosity value at the reference depth (:math:`\\eta_0`) in Pa.s.
    :param float zeta: Depth scale (:math:`y_{\\eta}`) in m.
    :param float zref: Reference depth (:math:`y_{\\text{ref}}`) in m.
    :param str model_name: Name of the model to which the viscosity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the viscosity is applied. Default is 0.

    .. note::
      The viscosity is given by:
        :math:`\\eta(y) = \\eta_0 \\exp \\left(- \\frac{(y_{\\text{ref}} - y)}{y_{\\eta}} \\right)`

    Attributes
    ----------

    .. py:attribute:: eta0
      :type: float

        Viscosity value at the reference depth (:math:`\\eta_0`) in :math:`\\text{Pa.s}`.

    .. py:attribute:: zeta
      :type: float

        Depth scale (:math:`y_{\\eta}`) in :math:`\\text{m}`.

    .. py:attribute:: zref
      :type: float

        Reference depth (:math:`y_{\\text{ref}}`) in :math:`\\text{m}`.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the viscosity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the viscosity is applied. Default is 0.
  """
  def __init__(self, eta0:float, zeta:float, zref:float, model_name:str="model_GENE3D", region:int=0) -> None:
    self.visc_type = 2
    self.eta0           = eta0
    self.zeta           = zeta
    self.zref           = zref
    Viscosity.__init__(self,model_name,region)

  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'eta(y) = eta_0 * exp( -(zref - y)/zeta )\n'
    s += f'\tRegion:   {self.region}\n'
    s += f'\tViscosity:{self.eta0}\n'
    s += f'\tZeta:     {self.zeta}\n'
    s += f'\tZref:     {self.zref}\n'
    return s 
  
class ViscosityArrhenius(Viscosity):
  """
  .. py:class:: ViscosityArrhenius(rock_name:str, Vmol:float=0.0, Tref:float=273.15, model_name:str="model_GENE3D", region:int=0, **kwargs)
  
    Class to apply an Arrhenius viscosity to a region of the model.
    For simplicity and to avoid potential errors, a database of standard rocks arrhenius flow laws is provided in
    the file arrhenius_flow_law.json.
    The file is loaded by the :class:`ViscosityArrhenius` class (and its children) 
    and the available rocks are stored in the class attribute :attr:`ViscosityArrhenius.__rocks__`.
    The dictionnary containing the rocks can be accessed with the class method :meth:`ViscosityArrhenius.arrhenius_flow_laws`.
  
    :param str rock_name: Name of the rock to retrieve the Arrhenius flow law parameters.
    :param float Vmol: Molar volume of the rock in  m\\ :sup:`3`.mol\\ :sup:`-1` . Default is 0.0.
    :param float Tref: Reference value to convert the temperature from Celsius to Kelvin. Default is 273.15.
    :param str model_name: Name of the model to which the viscosity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the viscosity is applied. Default is 0.
    :param kwargs: Additional arguments to erase parameters from the database if the rock is available or to provide the parameters if the rock is not in the database.
    
    .. note::
      The viscosity is given by:
        :math:`\\eta(\\mathbf u, p, T) = \\frac{1}{4} \\varepsilon^{(\\frac{1}{n} - 1)} \\left(\\frac{3}{4} A \\right)^{-\\frac{1}{n}} \\exp \\left( \\frac{E + p V}{n R T} \\right)`
      
      with :math:`\\varepsilon` the strain rate norm, 
      :math:`A` the preexpA parameter, 
      :math:`E` the entalpy, 
      :math:`n` the nexp parameter, 
      :math:`V` the molar volume, 
      :math:`R` the gas constant and 
      :math:`T` the temperature.
    
      
    Attributes
    ----------

    .. py:attribute:: preexpA
      :type: float

        Pre-exponential factor (:math:`A`) in Pa\\ :sup:`-n`.s\\ :sup:`-1`.

    .. py:attribute:: Ascale
      :type: float

        Scaling factor for the viscosity. 
        Should be :math:`10^6` if :math:`A` is given in MPa\\ :sup:`-n`.s\\ :sup:`-1` 
        and :math:`1` if :math:`A` is given in Pa\\ :sup:`-n`.s\\ :sup:`-1`.
      
    .. py:attribute:: entalpy
      :type: float

        Activation energy (:math:`E`) in J.mol\\ :sup:`-1`.

    .. py:attribute:: Vmol
      :type: float

        Molar volume of the rock in m\\ :sup:`3`.mol\\ :sup:`-1`.

    .. py:attribute:: nexp
      :type: float

        Exponent of the viscosity law (:math:`n`).

    .. py:attribute:: Tref
      :type: float

        Reference value to convert the temperature from Celsius to Kelvin. Default is 273.15.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the viscosity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the viscosity is applied. Default is 0.
  """
  # Define the standard rocks arrhenius flow laws as class attribute (not class instance attribute)
  # This ensures that the dictionnary is shared among all instances of the class and that if the dictionnary is modified, all instances will see the change
  with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"arrhenius_flow_laws.json"),"r") as f:
    __rocks__ = json.load(f)

  def __init__(self, rock_name:str, Vmol:float=0.0, Tref:float=273.15, model_name:str="model_GENE3D", region:int=0, **kwargs) -> None:
    self.visc_type = 3
    # generate the dictionnary of standard rocks arrhenius flow laws
    rock_param = ViscosityArrhenius.arrhenius_flow_laws()
    # verify if the requested rock is in the dictionnary
    if rock_name in rock_param:
      self.preexpA = kwargs.get("preexpA", rock_param[rock_name]["preexpA"])
      self.Ascale  = kwargs.get("Ascale",  rock_param[rock_name]["Ascale"])
      self.entalpy = kwargs.get("entalpy", rock_param[rock_name]["entalpy"])
      self.nexp    = kwargs.get("nexp",    rock_param[rock_name]["nexp"])
    else:
      if len(kwargs) == 0:
        s =  f"{self.__class__.__name__}: The requested flow law named: {rock_name} is not in the list of available flow laws\n"
        s += f"To work around, you must provide the following arguments to the constructor: preexpA, nexp, entalpy, Ascale.\n"
        s += f"Available flow laws are:\n"
        for r in rock_param:
          s += f"\t- {r}\n"
        raise ValueError(s)
      self.preexpA = kwargs["preexpA"]
      self.Ascale  = kwargs["Ascale"]
      self.entalpy = kwargs["entalpy"]
      self.nexp    = kwargs["nexp"]
    
    self.Vmol           = Vmol
    self.Tref           = Tref
    Viscosity.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'eta(u,p,T) = 0.25*Ascale*strainrate^(1/nexp - 1)*(0.75*preexpA)^(-1/nexp)*exp((entalpy + p*Vmol)/(nexp*R*T))\n'
    s += f'\tRegion:    {self.region}\n'
    s += f'\tpreexpA:   {self.preexpA}\n'
    s += f'\tAscale:    {self.Ascale}\n'
    s += f'\tentalpy:   {self.entalpy}\n'
    s += f'\tVmol:      {self.Vmol}\n'
    s += f'\tnexp:      {self.nexp}\n'
    s += f'\tTref:      {self.Tref}\n'
    return s
  
  @classmethod
  def arrhenius_flow_laws(cls) -> dict:
    """
    arrhenius_flow_laws(cls)
    Return the dictionnary of standard rocks arrhenius flow laws.
    The dictionnary is loaded from the file arrhenius_flow_laws.json.
    The file is loaded by the :class:`ViscosityArrhenius` class (and its children)
    and the available rocks are stored in the class attribute :attr:`ViscosityArrhenius.__rocks__`.
    
    :return: The dictionnary of standard rocks arrhenius flow laws.
    :rtype: dict
    """
    return cls.__rocks__

class ViscosityArrhenius2(ViscosityArrhenius):
  """
  .. py:class:: ViscosityArrhenius2(rock_name:str, Vmol:float=0.0, Tref:float=273.15, model_name:str="model_GENE3D", region:int=0, **kwargs)

    Class to apply an Arrhenius viscosity to a region of the model.
    
    .. note:: 
      This class is similar to the :class:`ViscosityArrhenius` class but the viscosity is given by:
        :math:`\\eta(\\mathbf u, p, T) = \\varepsilon^{(\\frac{1}{n} - 1)} A^{-\\frac{1}{n}} \\exp \\left( \\frac{E + p V}{n R T} \\right)`

  """
  def __init__(self, rock_name:str, Vmol:float=0.0, Tref:float=273.15, model_name:str="model_GENE3D", region:int=0, **kwargs) -> None:
    ViscosityArrhenius.__init__(self,rock_name,Vmol,Tref,model_name,region,**kwargs)
    self.visc_type = 4
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'eta(u,p,T) = Ascale*strainrate^(1/nexp - 1)*preexpA^(-1/nexp)*exp((entalpy + p*Vmol)/(nexp*R*T))\n'
    s += f'\tRegion:   {self.region}\n'
    s += f'\tpreexpA:  {self.preexpA}\n'
    s += f'\tAscale:   {self.Ascale}\n'
    s += f'\tentalpy:  {self.entalpy}\n'
    s += f'\tVmol:     {self.Vmol}\n'
    s += f'\tnexp:     {self.nexp}\n'
    s += f'\tTref:     {self.Tref}\n'
    return s

class ViscosityArrheniusDislDiff(ViscosityArrhenius):
  """
  .. py:class:: ViscosityArrheniusDislDiff(rock_name:str, preexpA_diff:float, Ascale_diff:float, entalpy_diff:float, Vmol_diff:float, pexp_diff:float, gsize:float, Vmol_disl:float=0.0, Tref:float=273.15, model_name:str="model_GENE3D", region:int=0, **kwargs)

    Class to apply a dislocation-diffusion viscosity to a region of the model.

    :param str rock_name: Name of the rock to retrieve the dislocation flow law parameters.
    :param float preexpA_diff: Pre-exponential factor for the diffusion flow law in Pa.s\\ :sup:`-1`.
    :param float Ascale_diff: Scaling factor for the diffusion viscosity. 
      Should be :math:`10^6` if :math:`A` is given in MPa.s\\ :sup:`-1`
      and :math:`1` if :math:`A` is given in Pa.s\\ :sup:`-1`.
    :param float entalpy_diff: Activation energy for the diffusion flow law in J.mol\\ :sup:`-1`.
    :param float Vmol_diff: Molar volume of the diffusion in m\\ :sup:`3`.mol\\ :sup:`-1`.
    :param float pexp_diff: Grain size exponent of the diffusion flow law.
    :param float gsize: Grain size in m.
    :param float Vmol_disl: Molar volume of the dislocation flow law in m\\ :sup:`3`.mol\\ :sup:`-1`. Default is 0.0.
    :param float Tref: Reference value to convert the temperature from Celsius to Kelvin. Default is 273.15.
    :param str model_name: Name of the model to which the viscosity is applied. Default is "model_GENE3D".
    :param int region: Region number to which the viscosity is applied. Default is 0.
    :param kwargs: Additional arguments to erase parameters from the database if the rock is available or to provide the parameters if the rock is not in the database.
    
    .. note::
      The viscosity is given by:

      .. math:: 
        \\eta_{\\text{disl}}(\\mathbf u, p, T) &= \\eta(\\mathbf u, p, T) = \\varepsilon^{(\\frac{1}{n} - 1)} A^{-\\frac{1}{n}} \\exp \\left( \\frac{E + p V}{n R T} \\right) \\\\
        \\eta_{\\text{diff}}(p,T) &= A \\exp \\left( \\frac{E + p V}{R T} \\right) g^{p} \\\\
        \\eta &= \\left( \\frac{1}{\\eta_{\\text{disl}} + \\frac{1}{\\eta_{\\text{diff}}}} \\right)^{-1}
    
    Attributes
    ----------

    .. py:attribute:: preexpA_diff
      :type: float
        
          Pre-exponential factor for the diffusion flow law in Pa.s\\ :sup:`-1`.

    .. py:attribute:: Ascale_diff
      :type: float

        Scaling factor for the diffusion viscosity. 
        Should be :math:`10^6` if :math:`A` is given in MPa.s\\ :sup:`-1`
        and :math:`1` if :math:`A` is given in Pa.s\\ :sup:`-1`.

    .. py:attribute:: entalpy_diff
      :type: float

        Activation energy for the diffusion flow law in J.mol\\ :sup:`-1`.

    .. py:attribute:: Vmol_diff

      Molar volume of the diffusion in m\\ :sup:`3`.mol\\ :sup:`-1`.

    .. py:attribute:: pexp_diff
      :type: float
        
          Grain size exponent of the diffusion flow law.

    .. py:attribute:: gsize
      :type: float

        Grain size in m.


    .. py:attribute:: preexpA_disl
      :type: float

        Pre-exponential factor for the dislocation flow law in Pa.s\\ :sup:`-1`.

    .. py:attribute:: Ascale_disl
      :type: float
        
          Scaling factor for the dislocation viscosity. 
          Should be :math:`10^6` if :math:`A` is given in MPa.s\\ :sup:`-1`
          and :math:`1` if :math:`A` is given in Pa.s\\ :sup:`-1`.

    .. py:attribute:: entalpy_disl
      :type: float

        Activation energy for the dislocation flow law in J.mol\\ :sup:`-1`.
      
    .. py:attribute:: Vmol_disl
      :type: float

        Molar volume of the dislocation in m\\ :sup:`3`.mol\\ :sup:`-1`.
        
    .. py:attribute:: nexp_disl
      :type: float

        Exponent of the dislocation viscosity law (:math:`n`).

    .. py:attribute:: Tref
      :type: float
        
          Reference value to convert the temperature from Celsius to Kelvin. Default is 273.15.

    .. py:attribute:: model_name
      :type: str

        Name of the model to which the viscosity is applied. Default is "model_GENE3D"

    .. py:attribute:: region
      :type: int

        Region number to which the viscosity is applied. Default is 0.
  """
  def __init__(self, rock_name:str,
               preexpA_diff:float, Ascale_diff:float, entalpy_diff:float, Vmol_diff:float, pexp_diff:float, gsize:float,
               Vmol_disl:float=0.0, Tref:float=273.15, model_name:str="model_GENE3D", region:int=0, **kwargs) -> None:
    ViscosityArrhenius.__init__(self,rock_name,Vmol_disl,Tref,model_name,region,**kwargs)
    # erase the attributes from the parent class
    self.preexpA_disl = self.preexpA
    self.Ascale_disl  = self.Ascale
    self.entalpy_disl = self.entalpy
    self.Vmol_disl    = self.Vmol
    self.nexp_disl    = self.nexp
    # delete the attributes from the parent class to avoid printing them in the sprint_option method
    del self.preexpA
    del self.Ascale
    del self.entalpy
    del self.Vmol
    del self.nexp

    self.preexpA_diff = preexpA_diff
    self.Ascale_diff  = Ascale_diff
    self.entalpy_diff = entalpy_diff
    self.Vmol_diff    = Vmol_diff
    self.pexp_diff    = pexp_diff
    self.gsize        = gsize
    self.visc_type    = 5
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}\n'
    s += f'eta_disl(u,p,T) = Ascale_disl*strainrate^(1/nexp_disl - 1)*preexpA_disl^(-1/nexp_disl)*exp((entalpy_disl + p*Vmol_disl)/(nexp_disl*R*T))\n'
    s += f'eta_diff(p,T)   = Ascale_diff*preexpA_diff^(-1)*gsize^pexp_diff*exp((entalpy_diff + p*Vmol_diff)/(R*T))\n'
    s += f'eta = (1/eta_disl + 1/eta_diff)^(-1)\n'
    s += f'\tRegion:               {self.region}\n'
    s += f'\tpreexpA dislocation:  {self.preexpA_disl}\n'
    s += f'\tAscale dislocation:   {self.Ascale_disl}\n'
    s += f'\tentalpy dislocation:  {self.entalpy_disl}\n'
    s += f'\tVmol dislocation:     {self.Vmol_disl}\n'
    s += f'\tnexp dislocation:     {self.nexp_disl}\n'
    s += f'\tpreexpA diffusion:    {self.preexpA_diff}\n'
    s += f'\tAscale diffusion:     {self.Ascale_diff}\n'
    s += f'\tentalpy diffusion:    {self.entalpy_diff}\n'
    s += f'\tVmol diffusion:       {self.Vmol_diff}\n'
    s += f'\tpexp diffusion:       {self.pexp_diff}\n'
    s += f'\tgsize:                {self.gsize}\n'
    s += f'\tTref:                 {self.Tref}\n'
    return s
