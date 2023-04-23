import sif_parser
import matplotlib.pyplot as plt
from scipy.special import voigt_profile as voigt
import numpy as np
import scipy.integrate as integrate
from functools import partial
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
    temperatures, left_peaks,left_peaks_err, right_peaks,right_peaks_err, left_peak_param, right_peak_param, fitting_function
):
    fig = plt.figure(figsize=(8,4))
    ax= fig.subplots(1,2)

    temp=np.linspace(0,260,200)
    ax[0].errorbar(temperatures,left_peaks,yerr=left_peaks_err, fmt='.k')
    ax[0].plot(temp,fitting_function(temp,*left_peak_param), 'r-', linewidth=2.0)
    ax[0].set_xlabel('temperature [K]')
    ax[0].set_ylabel('Peak Position [eV]')

    ax[1].errorbar(temperatures,right_peaks,yerr=right_peaks_err,fmt='.k')
    ax[1].plot(temp,fitting_function(temp,*right_peak_param), 'r-', linewidth=2.0)
    ax[1].set_xlabel('temperature [K]')
    ax[1].set_ylabel('Peak Position [eV]')

    fig.tight_layout()
    return fig

def plot_temperature_width_fig(
    temperatures, left_peak_widths,left_peak_widths_err, right_peak_widths, right_peak_widths_err, left_peak_width_param, right_peak_width_param, fitting_function
):
    fig = plt.figure(figsize=(8,4))
    ax= fig.subplots(1,2)

    temp=np.linspace(0,260,200)
    ax[0].errorbar(temperatures,left_peak_widths,yerr=left_peak_widths_err, fmt='.k')
    ax[0].plot(temp,fitting_function(temp,*left_peak_width_param), 'r-', linewidth=2.0)
    ax[0].set_xlabel('temperature [K]')
    ax[0].set_ylabel('Peak width [eV]')

    ax[1].errorbar(temperatures,right_peak_widths,yerr=right_peak_widths_err, fmt='.k')
    ax[1].plot(temp,fitting_function(temp,*right_peak_width_param), 'r-', linewidth=2.0)
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

def peak_temperature_function(
    T,E0,a,T_D=760
):
    quad=np.array(
        list(map(partial(integrate.quad,f1,0),T_D/T))
    )[:,0]
    return E0+a*((T/T_D)**4)*quad

def f1(x):
    return (x**3)/(np.exp(x)-1)

def f2(x):
    return (x**7)/((np.exp(x)-1)**2)

def peak_width_temperature_function(
    T,G0,a,T_D=760
):
    quad=np.array(
        list(map(partial(integrate.quad,f2,0),T_D/T))
    )[:,0]
    return G0+a*((T/T_D)**7)*quad
