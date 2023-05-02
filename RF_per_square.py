from RF_plots_season import read_files,path_chooser

emission_point = 1
f_string = path_chooser("summer",250)
global_flux = read_files(f_string, emission_point, RF=True, Atilla=False, O3=False)

print(global_flux[emission_point-1])