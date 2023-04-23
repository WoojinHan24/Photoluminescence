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
slit_widthes=[]
indexes=[]

for file in os.listdir(folder_name):
    if fnmatch.fnmatch(file, 'Ruby s?-?.sif'):
        data,wavelength = pl.get_sif_file(folder_name+file)
        photon_numbers.append(np.array(data))
        wavelength=wavelength*1e-9
        wavelengthes.append(np.array(wavelength))
        energy = h*c/wavelength/e
        energies.append(np.array(energy))
        slit_widthes.append(int(file[6]))
        indexes.append(int(file[-5]))

    if fnmatch.fnmatch(file, 'Ruby s??-?.sif'):
        data,wavelength = pl.get_sif_file(folder_name+file)
        photon_numbers.append(np.array(data))
        wavelength=wavelength*1e-9
        wavelengthes.append(np.array(wavelength))
        energy = h*c/wavelength/e
        energies.append(np.array(energy))
        slit_widthes.append(int(file[6:8]))
        indexes.append(int(file[-5]))
    
for photon_number,wavelength,energy,slit_width,index in zip(photon_numbers,wavelengthes,energies,slit_widthes,indexes): 
    raw_fig=pl.plot_raw_fig(photon_number,wavelength)
    raw_fig.savefig(f"./results/Ruby s{slit_width}-{index}_raw_fig.png")


