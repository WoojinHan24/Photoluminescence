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
e = constants.e
photon_numbers=[]
wavelengthes=[]
energies=[]
temperatures=[]
left_peaks=[]
right_peaks=[]
left_peak_widths_g=[]
right_peak_widths_g=[]
optical_broadening=[]
left_peak_heights=[]
right_peak_heights=[]
parameter_covariances=[]

for file in os.listdir(folder_name):
    if fnmatch.fnmatch(file, 'Rby *.sif'):
        data,wavelength = pl.get_sif_file(folder_name+file)
        photon_numbers.append(np.array(data))
        temperatures.append(np.float64(file[4:-4]))
        wavelength=wavelength*1e-9
        wavelengthes.append(np.array(wavelength))
        energy = h*c/wavelength/e
        energies.append(np.array(energy))


for photon_number,temperature,wavelength,energy in zip(photon_numbers,temperatures,wavelengthes,energies): 

    raw_fig=pl.plot_raw_fig(photon_number,wavelength)
    raw_fig.savefig(f"./results/Ruby({temperature})_raw_fig.png")

    if temperature<30:
        #low temperature one peak case walking
        p1=np.array([1.783,1e-3,1e-3,3,4400])
        param, param_covariance=opt.curve_fit(
            pl.one_voigt_function,
            energy.flatten(),
            photon_number.flatten(),
            p1,
            maxfev=60000
        )
        fitted_fig=pl.plot_param_fig(photon_number,energy,param,pl.one_voigt_function_param)
        fitted_fig.savefig(f"./results/Ruby({temperature})_voigt_fitted_fig.png")
        left_peaks.append(param[0])
        right_peaks.append(np.nan)
        optical_broadening.append(param[1])
        left_peak_widths_g.append(param[2])
        right_peak_widths_g.append(np.nan)
        left_peak_heights.append(param[3])
        right_peak_heights.append(np.nan)
        parameter_covariances.append(param_covariance)

        continue

    p0=np.array([1.783,1.790,1e-3,1e-4,1e-4,3,1,4400])
    param, param_covariance=opt.curve_fit(
        pl.two_voigt_function,
        energy.flatten(),
        photon_number.flatten(),
        p0,
        maxfev=60000
    )

    fitted_fig=pl.plot_param_fig(photon_number,energy,param,pl.two_voigt_function_param)
    fitted_fig.savefig(f"./results/Ruby({temperature})_voigt_fitted_fig.png")
    left_peaks.append(param[0])
    right_peaks.append(param[1])
    optical_broadening.append(param[2])
    left_peak_widths_g.append(param[3])
    right_peak_widths_g.append(param[4])
    left_peak_heights.append(param[5])
    right_peak_heights.append(param[6])
    parameter_covariances.append(param_covariance)
    



temperature_two_peak_fig=pl.plot_temperature_peak_fig(temperatures,left_peaks,right_peaks)
temperature_width_fig=pl.plot_temperature_width_fig(temperatures,left_peak_widths_g,right_peak_widths_g)

temperature_two_peak_fig.savefig(f"./results/Ruby_temperature_peak_fig.png")
temperature_width_fig.savefig(f"./results/Ruby_temperature_width_fig.png")


