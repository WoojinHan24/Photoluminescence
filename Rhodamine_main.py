import photoluminescence as pl
import fnmatch
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants
import scipy.optimize as opt


folder_name="./raw_data/"

h = constants.Planck
c = constants.c
e = constants.e
photon_numbers=[]
wavelengthes=[]
energies=[]
indexes=[]


for file in os.listdir(folder_name):
    if fnmatch.fnmatch(file, 'Rt s25 05-20 *.sif'):
        data,wavelength = pl.get_sif_file(folder_name+file)
        photon_numbers.append(np.array(data))
        wavelength=wavelength*1e-9
        wavelengthes.append(np.array(wavelength))
        energy = h*c/wavelength/e
        energies.append(np.array(energy))
        indexes.append(int(file[-5]))
    
total_fig = plt.figure()
ax = total_fig.add_subplot(1,1,1)
color = [' ','k-', 'r-', 'b-', 'g-','y-']
for photon_number,wavelength,energy,index in zip(photon_numbers,wavelengthes,energies,indexes): 
    raw_fig=pl.plot_raw_fig(photon_number,wavelength)
    raw_fig.savefig(f"./results/R6G({index})_raw_fig.png")

    p3=np.array([2.23,0.1,7000,10000])
    param, param_covariance=opt.curve_fit(
        pl.one_gaussian_function,
        energy.flatten(),
        photon_number.flatten(),
        p3,
        maxfev=600000000
    )
    print(index,param)
    fitted_fig=pl.plot_param_fig(photon_number,energy,param,pl.one_gaussian_function)
    fitted_fig.savefig(f"./results/R6G({index})_gaussian_fitted_fig.png")

    ax.plot(energy,photon_number,color[index])


ax.set_xlabel("photon energy [$eV$]")
ax.set_ylabel("detected photon number")
total_fig.savefig("./results/R6G_total_fig.png")