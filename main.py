import photoluminescence as pl
import fnmatch
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants
import scipy.optimize as opt


folder_name="./raw_data_temperature/"

h = constants.Planck
c = constants.c
photon_numbers=[]
wavelengthes=[]
energies=[]
temperatures=[]
for file in os.listdir(folder_name):
    if fnmatch.fnmatch(file, 'Rby *.sif'):
        data,wavelength = pl.get_sif_file(folder_name+file)
        photon_numbers.append(np.array(data))
        temperatures.append(np.float64(file[4:-4]))
        wavelength=wavelength*1e-9
        wavelengthes.append(np.array(wavelength))
        energy = h*c/wavelength
        energies.append(np.array(energy))


for photon_number,temperature,wavelength,energy in zip(photon_numbers,temperatures,wavelengthes,energies):
    print(type(energy))
    print(type(photon_number))   
    #param=opt.curve_fit(pl.three_Lorentzian_function,energy.flatten(),photon_number.flatten(),np.array([3*1e-19,4*1e-19,2*1e-19,1e-20,1e-20,1e-20,1400,1600,400,4200]))
    param, param_covariance=opt.curve_fit(
        pl.three_Lorentzian_function,
        energy.flatten(),
        photon_number.flatten(),
        np.array([2.86*1e-19,2.869*1e-19,2.95*1e-19,1e-22,1e-22,5e-21,1500,1000,100,4400]),
        maxfev=6000
    )
    
    raw_fig=pl.plot_raw_fig(photon_number,energy,param)
    raw_fig.savefig(f"./results/Ruby({temperature})_rawfig.png")
    print(param)

