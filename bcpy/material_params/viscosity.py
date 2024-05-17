from bcpy import MaterialConstants

class Viscosity(MaterialConstants):
  def __init__(self, model_name:str, region:int) -> None:
    MaterialConstants.__init__(self,model_name,region)

class ViscosityConstant(Viscosity):
  def __init__(self, model_name:str, region:int, viscosity:float) -> None:
    self.viscosity_type = 0
    self.eta0           = viscosity
    Viscosity.__init__(self,model_name,region)
  
  def __str__(self) -> str:
    s = f'{self.__class__.__name__}:\n'
    s += f'\tRegion:    {self.region}\n'
    s += f'\tViscosity: {self.eta0}\n'
    return s
  
class ViscosityFrankK(Viscosity):
  def __init__(self, model_name: str, region: int, eta0:float, exponent:float) -> None:
    self.viscosity_type = 1
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
  def __init__(self, model_name: str, region: int, eta0:float, zeta:float, zref:float) -> None:
    self.viscosity_type = 2
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
  def __init__(self, model_name: str, region: int, rock_name:str, Vmol:float=0.0, Tref:float=273.15, **kwargs) -> None:
    self.viscosity_type = 3
    # generate the dictionnary of standard rocks arrhenius flow laws
    rock_param = self.arrhenius_flow_laws()
    # verify if the requested rock is in the dictionnary
    if rock_name in rock_param:
      self.preexpA = rock_param[rock_name]["preexpA"]
      self.Ascale  = rock_param[rock_name]["Ascale"]
      self.entalpy = rock_param[rock_name]["entalpy"]
      self.nexp    = rock_param[rock_name]["nexp"]
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
    s += f'\tRegion:   {self.region}\n'
    s += f'\tpreexpA:  {self.preexpA}\n'
    s += f'\tAscale:   {self.Ascale}\n'
    s += f'\tentalpy:  {self.entalpy}\n'
    s += f'\tVmol:     {self.Vmol}\n'
    s += f'\tnexp:     {self.nexp}\n'
    s += f'\tTref:     {self.Tref}\n'
    return s
  
  def arrhenius_flow_laws(self) -> None:
    rock_param = dict()
    rock_param["Quartz"]            = {"preexpA":1.0e-3, "nexp":2.0, "entalpy":167.0e3, "Ascale":1.0e6}
    rock_param["Plagioclase(An75)"] = {"preexpA":3.3e-4, "nexp":3.2, "entalpy":238.0e3, "Ascale":1.0e6}
    rock_param["Orthopyroxene"]     = {"preexpA":3.2e-1, "nexp":2.4, "entalpy":293.0e3, "Ascale":1.0e6}
    rock_param["Clinopyroxene"]     = {"preexpA":15.7,   "nexp":2.6, "entalpy":335.0e3, "Ascale":1.0e6}
    rock_param["Granite"]           = {"preexpA":1.8e-9, "nexp":3.2, "entalpy":123.0e3, "Ascale":1.0e6}
    rock_param["Granite(wet)"]      = {"preexpA":2.0e-4, "nexp":1.9, "entalpy":137.0e3, "Ascale":1.0e6}
    rock_param["Quartzite"]         = {"preexpA":6.7e-6, "nexp":2.4, "entalpy":156.0e3, "Ascale":1.0e6}
    rock_param["Quartzite(wet)"]    = {"preexpA":3.2e-4, "nexp":2.3, "entalpy":154.0e3, "Ascale":1.0e6}
    rock_param["QuartzDiorite"]     = {"preexpA":1.3e-3, "nexp":2.4, "entalpy":219.0e3, "Ascale":1.0e6}
    rock_param["Diabase"]           = {"preexpA":2.0e-4, "nexp":3.4, "entalpy":260.0e3, "Ascale":1.0e6}
    rock_param["Anorthosite"]       = {"preexpA":3.2e-4, "nexp":3.2, "entalpy":238.0e3, "Ascale":1.0e6}
    rock_param["FelsicGranulite"]   = {"preexpA":8.0e-3, "nexp":3.1, "entalpy":243.0e3, "Ascale":1.0e6}
    rock_param["MaficGranulite"]    = {"preexpA":1.4e4,  "nexp":4.2, "entalpy":445.0e3, "Ascale":1.0e6}
    rock_param["Peridotite(dry)"]   = {"preexpA":2.5e4,  "nexp":3.5, "entalpy":532.0e3, "Ascale":1.0e6}
    rock_param["Peridotite(wet)"]   = {"preexpA":2.0e3,  "nexp":4.0, "entalpy":471.0e3, "Ascale":1.0e6}
    return rock_param

class ViscosityArrhenius2(ViscosityArrhenius):
  def __init__(self, model_name: str, region: int, rock_name:str, Vmol:float=0.0, Tref:float=273.15, **kwargs) -> None:
    ViscosityArrhenius.__init__(self,model_name,region,rock_name,Vmol,Tref,**kwargs)
    self.viscosity_type = 4
  
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
  def __init__(self, model_name: str, region: int, rock_name:str,
               preexpA_diff:float, Ascale_diff:float, entalpy_diff:float, Vmol_diff:float, pexp_diff:float, gsize:float,
               Vmol_disl:float=0.0, Tref:float=273.15, **kwargs) -> None:
    ViscosityArrhenius.__init__(self,model_name,region,rock_name,Vmol_disl,Tref,**kwargs)
    # erase the attributes from the parent class
    self.preexpA_disl = self.preexpA
    self.Ascale_disl = self.Ascale
    self.entalpy_disl = self.entalpy
    self.Vmol_disl = self.Vmol
    self.nexp_disl = self.nexp
    # delete the attributes from the parent class to avoid printing them in the sprint_option method
    del self.preexpA
    del self.Ascale
    del self.entalpy
    del self.Vmol
    del self.nexp

    self.preexpA_diff   = preexpA_diff
    self.Ascale_diff    = Ascale_diff
    self.entalpy_diff   = entalpy_diff
    self.Vmol_diff      = Vmol_diff
    self.pexp_diff      = pexp_diff
    self.gsize          = gsize
    self.viscosity_type = 5
  
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