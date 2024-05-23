class MaterialConstants:
  def __init__(self, model_name:str="model_GENE3D", region:int=0) -> None:
    self.model_name = model_name
    self.region     = region
  
  def sprint_option(self):
    attributes = vars(self)
    s = ""
    for p in attributes:
      if p in ['model_name','region']:
        continue
      if type(attributes[p]) is str: fmt = f"{attributes[p]}"
      else:                          fmt = f"{attributes[p]:g}"
      s += f"-{self.model_name}_{p}_{self.region} {fmt}\n"
    return s