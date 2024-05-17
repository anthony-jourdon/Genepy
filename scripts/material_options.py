import numpy as np
import bcpy as bp

def main():
  region = 1
  region_0 = []
  region_0.append(bp.DensityBoussinesq(2900.0,3e-5,1e-11,"model_GENE3D",region))
  #region_0.append(bp.DensityTable(2700.0,"density_map.txt","model_GENE3D",region))
  region_0.append(bp.PlasticDruckerPrager(np.deg2rad(30.0),np.deg2rad(5.0),2.0e7,5.0e6,1.0e6,4.0e8,"model_GENE3D",region))
  region_0.append(bp.SofteningLinear(0.0,0.5,"model_GENE3D",region))
  region_0.append(bp.ViscosityArrhenius2("Granite",model_name="model_GENE3D",region=region))
  #region_0.append(bp.ViscosityArrhenius2("Salt",model_name="model_GENE3D",region=region,preexpA=3.0,nexp=2.1,entalpy=25.0e3,Ascale=1.0e6))
  #region_0.append(bp.ViscosityArrheniusDislDiff("Granite",model_name="model_GENE3D",region=region,preexpA_diff=1.0e-1,Ascale_diff=1.0e6,entalpy_diff=25.0e3,Vmol_diff=0.0,pexp_diff=2.5,gsize=25.0e3))

  opt = "Material options:\n"
  for r in region_0:
    opt += r.sprint_option()
  print(opt)

  

if __name__ == "__main__":
  main()