import photoluminescence as pl
import fnmatch
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants
import scipy.optimize as opt
import warnings


warnings.filterwarnings(action='ignore')

folder_name="./raw_data_temperature/"

h = constants.Planck
c = constants.c
e = constants.e
photon_numbers=[]
wavelengthes=[]
energies=[]
temperatures=[]
left_peaks=[]
left_peaks_err=[]
right_peaks=[]
right_peaks_err=[]
left_peak_widths_g=[]
left_peak_widths_g_err=[]
right_peak_widths_g=[]
right_peak_widths_g_err=[]
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
        p1=np.array([1.783,5.7*1e-4,1.5*1e-5,2.2,4400])
        param, param_covariance=opt.curve_fit(
            pl.one_voigt_function,
            energy.flatten(),
            photon_number.flatten(),
            p1,
            maxfev=60000
        )
        fitted_fig=pl.plot_param_fig(photon_number,energy,param,pl.one_voigt_function)
        fitted_fig.savefig(f"./results/Ruby({temperature})_voigt_fitted_fig.png")
        left_peaks.append(param[0])
        left_peaks_err.append(param_covariance[0][0]+0.0001)
        right_peaks.append(np.nan)
        right_peaks_err.append(np.nan)
        optical_broadening.append(param[1])
        left_peak_widths_g.append(np.abs(param[2]))
        left_peak_widths_g_err.append(param_covariance[2][2]+0.0001)
        right_peak_widths_g.append(np.nan)
        right_peak_widths_g_err.append(np.nan)
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

    fitted_fig=pl.plot_param_fig(photon_number,energy,param,pl.two_voigt_function)
    fitted_fig.savefig(f"./results/Ruby({temperature})_voigt_fitted_fig.png")
    left_peaks.append(param[0])
    left_peaks_err.append(param_covariance[0][0]+0.0001)
    right_peaks.append(param[1])
    right_peaks_err.append(param_covariance[1][1]+0.0001)
    optical_broadening.append(param[2])
    left_peak_widths_g.append(np.abs(param[3]))
    left_peak_widths_g_err.append(param_covariance[3][3]+0.0001)
    right_peak_widths_g.append(np.abs(param[4]))
    right_peak_widths_g_err.append(param_covariance[4][4]+0.0001)
    left_peak_heights.append(param[5])
    right_peak_heights.append(param[6])
    parameter_covariances.append(param_covariance)



left_peak_widths=[np.sqrt(x**2+y**2) for x,y in zip(left_peak_widths_g,optical_broadening)]
right_peak_widths=[np.sqrt(x**2+y**2) for x,y in zip(right_peak_widths_g,optical_broadening)]

p_R1=np.array([1.78,-0.04])
left_peak_param, left_peak_param_cov = opt.curve_fit(
    pl.peak_temperature_function,
    temperatures,
    left_peaks,
    p_R1,
    maxfev=60000
)

p_R2=np.array([1.81,-0.04])
right_peak_param, right_peak_param_cov = opt.curve_fit(
    pl.peak_temperature_function,
    [temperature for temperature,right_peak in zip(temperatures,right_peaks) if np.isnan(right_peak)==False],
    [right_peak for right_peak in right_peaks if np.isnan(right_peak)==False],
    p_R2,
    maxfev=60000
)

p_R1_g=np.array([5.9*1e-7,2.31*1e-3])
left_peak_width_param, _ = opt.curve_fit(
    pl.peak_width_temperature_function,
    temperatures,
    left_peak_widths,
    p_R1_g,
    maxfev=60000
)

p_R2_g=np.array([0,2.31*1e-3])
right_peak_width_param, _ = opt.curve_fit(
    pl.peak_width_temperature_function,
    [temperature for temperature,right_peak_width in zip(temperatures,right_peak_widths) if np.isnan(right_peak_width)==False],
    [right_peak_width for right_peak_width in right_peak_widths if np.isnan(right_peak_width)==False],
    p_R2_g,
    maxfev=60000
)

print(left_peak_param)
print(right_peak_param)
print(left_peak_width_param)
print(right_peak_width_param)

temperature_two_peak_fig=pl.plot_temperature_peak_fig(temperatures,left_peaks,left_peaks_err,right_peaks,right_peaks_err,left_peak_param,right_peak_param,pl.peak_temperature_function)
temperature_width_fig=pl.plot_temperature_width_fig(temperatures,left_peak_widths,left_peak_widths_g_err,right_peak_widths,right_peak_widths_g_err,left_peak_width_param,right_peak_width_param,pl.peak_width_temperature_function)

temperature_two_peak_fig.savefig(f"./results/Ruby_temperature_peak_fig.png")
temperature_width_fig.savefig(f"./results/Ruby_temperature_width_fig.png")


R1=np.array([pl.one_voigt_function(left_peak,left_peak, S, left_peak_width, left_peak_height,0) for left_peak, S, left_peak_width, left_peak_height in zip(left_peaks, optical_broadening, left_peak_widths, left_peak_heights)])
R2=np.array([pl.one_voigt_function(right_peak,right_peak, S, right_peak_width, right_peak_height,0) for right_peak,right_peak, S, right_peak_width, right_peak_height in zip(right_peaks,right_peaks, optical_broadening, right_peak_widths, right_peak_heights)])
print(R1)
height_ratios=R2/R1
cnt_err=10
height_ratios_err = np.sqrt((1/R1)**2 + (R2/R1**2)**2)*cnt_err

temperature_height_ratio_fig=pl.plot_temperature_height_ratio_raw_plot(temperatures, height_ratios, height_ratios_err)
temperature_height_ratio_fig.savefig("./results/Ruby_temperature_height_ratio_fig.png")
temperature_height_ratio_loglog_fig=pl.plot_temperature_height_ratio_loglog_plot(np.array(temperatures),height_ratios,height_ratios_err)
temperature_height_ratio_loglog_fig.savefig("./results/Ruby_temperature_height_ratio_loglog_fig.png")