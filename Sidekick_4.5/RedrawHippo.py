import HippoAuto,os
wd = os.getcwd()
mutant = wd.split('/')[-2]
HippoAuto.plot_matrix(mutant,"energy_surface.png")
