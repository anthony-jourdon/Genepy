import numpy as np
import genepy as gp

def main():
  region = 1

  sources = (
    gp.EnergySourceConstant(1.5e-6,"model_GENE3D",region),
    gp.EnergySourceShearHeating("model_GENE3D",region)
  )

  energy_sources = gp.EnergySource(*sources,model_name="model_GENE3D",region=region)

  region_0 = [
    gp.DensityBoussinesq(2900.0,3e-5,1e-11,"model_GENE3D",region),
    gp.PlasticDruckerPrager(np.deg2rad(30.0),np.deg2rad(5.0),2.0e7,5.0e6,1.0e6,4.0e8,"model_GENE3D",region),
    gp.SofteningLinear(0.0,0.5,"model_GENE3D",region),
    gp.ViscosityArrhenius2("Granite",model_name="model_GENE3D",region=region),
    gp.Energy(heat_source=energy_sources,conductivity=2.9,model_name="model_GENE3D",region=region)
  ]

  opt = "Material options:\n"
  for r in region_0:
    opt += r.sprint_option()
  print(opt)

  

if __name__ == "__main__":
  main()