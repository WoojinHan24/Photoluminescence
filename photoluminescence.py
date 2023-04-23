import sif_parser
import matplotlib.pyplot as plt
from scipy.special import voigt_profile as voigt
import numpy as np
from scipy.constants import pi

def get_sif_file(
    file_name
):
    data,info = sif_parser.np_open(file_name)
    if data.shape != (1,1,1024) :
        print("data shape irregular ", data.shape, "is not (1,1,1024)")
        
    data=data[0][0]
    wavelength=sif_parser.utils.extract_calibration(info)
    return data,wavelength

def plot_raw_fig(
    photon_number,wavelength
):
    raw_fig = plt.figure(figsize=(6,4))
    ax=raw_fig.add_subplot(1,1,1)

    ax.plot(wavelength,photon_number,'k-', linewidth=1.0)
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('detected photon number')
    return raw_fig

def plot_param_fig(
    photon_number,energy,param, fitting_function
):
    param_fig = plt.figure(figsize=(6,4))
    ax=param_fig.add_subplot(1,1,1)

    ax.plot(energy,photon_number,'k-', linewidth=1.0)
    ax.set_xlabel('photon energy [eV]')
    ax.set_ylabel('detected photon number')

    ax.plot(energy,fitting_function(energy,*param),'r-')
    return param_fig

def plot_temperature_peak_fig(
    temperatures, left_peaks, right_peaks
):
    fig = plt.figure(figsize=(8,4))
    ax= fig.subplots(1,2)

    ax[0].plot(temperatures,left_peaks, 'ko')
    ax[0].set_xlabel('temperature [K]')
    ax[0].set_ylabel('Peak Position [eV]')

    ax[1].plot(temperatures,right_peaks, 'ko')
    ax[1].set_xlabel('temperature [K]')
    ax[1].set_ylabel('Peak Position [eV]')

    fig.tight_layout()
    return fig

def plot_temperature_width_fig(
    temperatures, left_peak_widths, right_peak_widths
):
    fig = plt.figure(figsize=(8,4))
    ax= fig.subplots(1,2)

    ax[0].plot(temperatures,left_peak_widths, 'ko')
    ax[0].set_xlabel('temperature [K]')
    ax[0].set_ylabel('Peak width [eV]')

    ax[1].plot(temperatures,right_peak_widths, 'ko')
    ax[1].set_xlabel('temperature [K]')
    ax[1].set_ylabel('Peak width [eV]')

    fig.tight_layout()
    return fig

def two_voigt_function(
    E, E1, E2, S, G1, G2, A1, A2, C
):
    #E is energy peak position,
    #S is the Gaussian width
    #G is the Lorentzian width
    #A is the amplitude of each peak
    #C is the constant counts
    return A1*voigt(E-E1,S,G1)+A2*voigt(E-E2,S,G2)+C

def one_voigt_function(
    E, E1, S, G1, A1, C
):
    #E is energy peak position,
    #S is the Gaussian width
    #G is the Lorentzian width
    #A is the amplitude of each peak
    #C is the constant counts
    return A1*voigt(E-E1,S,G1)+C

def one_gaussian_function(
    E,E0,S,A,C
):
    return A/np.sqrt(2*pi)/S * np.exp(-0.5*(E-E0)**2/S**2) +C