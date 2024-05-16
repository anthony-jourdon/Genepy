import numpy as np
import bcpy as bp

def main():
  region = 0
  region_0 = []
  region_0.append(bp.DensityBoussinesq("model_GENE3D",region,2900.0,3e-5,1e-11))
  #region_0.append(bp.DensityTable("model_GENE3D",region,2700.0,"density_map.txt"))
  region_0.append(bp.PlasticDruckerPrager("model_GENE3D",region,np.deg2rad(30.0),np.deg2rad(5.0),2.0e7,5.0e6,1.0e6,4.0e8))
  region_0.append(bp.SofteningLinear("model_GENE3D",region,0.0,0.5))
  region_0.append(bp.ViscosityArrhenius2("model_GENE3D",region,6.7e-6,1e6,156e3,0,2.4))

  opt = "Material options:\n"
  for r in region_0:
    opt += r.sprint_option()
  print(opt)

if __name__ == "__main__":
  main()