from genepy import Domain
import numpy as np

class GeoPoints:
  def __init__(self,mesh_size:float) -> None:
    self.points    = []
    self.mesh_size = mesh_size
    return

class GeoLayerBox(GeoPoints):
  def __init__(self,x_layers:np.ndarray,y_layers:np.ndarray,z_layers:np.ndarray,mesh_size:float) -> None:
    GeoPoints.__init__(self,mesh_size)
    self.x = x_layers
    self.y = y_layers
    self.z = z_layers
    self.n = np.array([self.x.shape[0], self.y.shape[0], self.z.shape[0]], dtype=np.int32)
    self.coordinates()
    return

  def coordinates(self):
    self.points.append(f"cl__1 = {self.mesh_size};")
    for k in range(self.n[2]):
      for j in range(self.n[1]):
        for i in range(self.n[0]):
          pidx = i + j*self.n[0] + k*self.n[0]*self.n[1]
          self.points.append(f'Point({pidx}) = '+'{'+ f'{self.x[i]}, {self.y[j]}, {self.z[k]}, cl__1' +'};')
    return
  
  def write(self,fname:str,delimiter:str='\n'):
    print(delimiter.join(self.points))
    with open(fname,'w') as f:
      f.write(delimiter.join(self.points))
    return

def test_box():
  Ox = 0.0
  Lx = 600.0e3
  nx = 2

  Oz = 0.0
  Lz = 300.0e3
  nz = 2

  Oy     = -250.0e3
  y_asth = -120.0e3
  y_lith = -35.0e3
  y_lc   = -20.0e3
  Ly     = 20.0e3

  x = np.linspace(Ox,Lx,nx)
  y = np.array([Oy, y_asth, y_lith, y_lc, Ly], dtype=np.float32)
  z = np.linspace(Oz,Lz,nz)

  lc = 2.0e+5
  box = GeoLayerBox(x,y,z,lc)
  box.write('box.geo')
  return

def test_subduction():
  return

if __name__ == "__main__":
  test_box()