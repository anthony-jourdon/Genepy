import numpy as np
from bcpy import domain

class MeshRefinement(domain.Domain):
  def __init__(self, Domain, refinement_params) -> None:
    domain.Domain.__init__(self,Domain.dim,Domain.O,Domain.L,Domain.n)
    self.params = refinement_params
    for d in self.params:
      if d not in ["x","y","z"]:
        raise ValueError(f"Invalid direction, possible values are \"x\", \"y\", \"z\", found: {d}")
      for p in self.params[d]:
        if p not in ["x_initial","x_refined"]:
          raise ValueError(f"Invalid parameter, possible values are \"x_initial\", \"x_refined\", found: {p}")

      self.params[d]["x_initial"] = np.asarray(self.params[d]["x_initial"],dtype=np.float64)
      self.params[d]["x_refined"] = np.asarray(self.params[d]["x_refined"],dtype=np.float64)
      if self.params[d]["x_initial"].shape != self.params[d]["x_refined"].shape:
        raise ValueError(f"Direction {d}: array in \"x_initial\" of shape {self.params[d]['x_initial'].shape} and \"x_refined\" {self.params[d]['x_refined'].shape} must have the same shape")
    return

  def __str__(self) -> str:
    s  = f"{self.__class__.__name__}:\n"
    for d in self.params:
      if   d == "x": dim = 0
      elif d == "y": dim = 1
      elif d == "z": dim = 2

      s += f"\tDirection: {d}\n"
      for p in self.params[d]:
        s += f"\t\t{p}           : {self.params[d][p]}\n"
        if p == "x_initial" or p == "x_refined":
          s += f"\t\t{p} normalized: {self.normalize(self.params[d][p],dim)}\n"
    return s
  
  def normalize(self,x,dim):
    return (x - self.O[dim]) / (self.L[dim] - self.O[dim])
  
  def refine_direction(self,dim,x_initial,x_refined):
    # convert tuple to list (tuples are immutable)
    coor = list(self.num_coor)
    # get all points in the dim direction
    xc_ini = coor[dim]
    # normalize
    xc_ini = self.normalize(xc_ini,dim)
    # linear interpolation
    xc_refined = np.interp(xc_ini,x_initial,x_refined)
    # rescale
    coor[dim] = self.O[dim] + xc_refined * (self.L[dim] - self.O[dim])
    # restore as a tuple
    self.num_coor = tuple(coor)
    return
  
  def refine(self):
    for d in self.params:
      if   d == "x": dim = 0
      elif d == "y": dim = 1
      elif d == "z": dim = 2

      x_initial = self.normalize(self.params[d]["x_initial"],dim)
      x_refined = self.normalize(self.params[d]["x_refined"],dim)
      self.refine_direction(dim,x_initial,x_refined)
    return

def test():
  from bcpy import writers
  O = np.array([0,-250,0],dtype=np.float64)
  L = np.array([600,0,600],dtype=np.float64)
  n = np.array([45,18,33],dtype=np.int32)
  d = domain.Domain(3,O,L,n)
  #w = writers.WriteVTS(d,vtk_fname="mesh_refinement.vts")

  refinement = {"x": {"x_initial": np.array([0,12,60,540,588,600], dtype=np.float64),
                      "x_refined": np.array([0,100,150,450,500,600], dtype=np.float64)},
                "y": {"x_initial": np.array([-250,-180,-87.5,0], dtype=np.float64),
                      "x_refined": np.array([-250,-50,-16.25,0], dtype=np.float64)},
                "z": {"x_initial": np.array([0,12,60,540,588,600], dtype=np.float64),
                      "x_refined": np.array([0,180,240,360,420,600], dtype=np.float64)}}

  m = MeshRefinement(d,refinement)
  print(m)
  m.refine()
  w = writers.WriteVTS(m,vtk_fname="mesh_refinement.vts")
  w.write_vts()

if __name__ == "__main__":
  test()