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
  def __init__(self, model_name: str, region: int, preexpA:float, Ascale:float, entalpy:float, Vmol:float, nexp:float, Tref:float=273.15) -> None:
    self.viscosity_type = 3
    self.preexpA        = preexpA
    self.Ascale         = Ascale
    self.entalpy        = entalpy
    self.Vmol           = Vmol
    self.nexp           = nexp
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

class ViscosityArrhenius2(Viscosity):
  def __init__(self, model_name: str, region: int, preexpA:float, Ascale:float, entalpy:float, Vmol:float, nexp:float, Tref:float=273.15) -> None:
    self.viscosity_type = 4
    self.preexpA        = preexpA
    self.Ascale         = Ascale
    self.entalpy        = entalpy
    self.Vmol           = Vmol
    self.nexp           = nexp
    self.Tref           = Tref
    Viscosity.__init__(self,model_name,region)
  
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

class ViscosityArrheniusDislDiff(Viscosity):
  def __init__(self, model_name: str, region: int, 
               preexpA_disl:float, Ascale_disl:float, entalpy_disl:float, Vmol_disl:float, nexp_disl:float,
               preexpA_diff:float, Ascale_diff:float, entalpy_diff:float, Vmol_diff:float, pexp_diff:float, gsize:float,
               Tref:float=273.15) -> None:
    self.viscosity_type = 5
    self.preexpA_disl   = preexpA_disl
    self.Ascale_disl    = Ascale_disl
    self.entalpy_disl   = entalpy_disl
    self.Vmol_disl      = Vmol_disl
    self.nexp_disl      = nexp_disl
    self.preexpA_diff   = preexpA_diff
    self.Ascale_diff    = Ascale_diff
    self.entalpy_diff   = entalpy_diff
    self.Vmol_diff      = Vmol_diff
    self.pexp_diff      = pexp_diff
    self.gsize          = gsize
    self.Tref           = Tref
    Viscosity.__init__(self,model_name,region)
  
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
    s += f'\tTref:     {self.Tref}\n'
    return s