import sif_parser
import matplotlib.pyplot as plt
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
    photon_number,energy,param
):
    raw_fig = plt.figure(figsize=(6,4))
    ax=raw_fig.add_subplot(1,1,1)

    ax.plot(energy,photon_number,'k-', linewidth=1.0)
    ax.set_xlabel('photon energy [J]')
    ax.set_ylabel('detected photon number')
    ax.plot(energy,three_Lorentzian_function_override(energy,param),'r-')
    return raw_fig

def three_Lorentzian_function(
    E,E1,E2,E3,D1,D2,D3,A1,A2,A3,C
):
    #E is energy peak position,
    #D is the width of each peak
    #A is amplitude of each lorentzian peak

    return Lorentzian_function(E,E1,D1,A1)+Lorentzian_function(E,E2,D2,A2)+Lorentzian_function(E,E3,D3,A3)+C


def three_Lorentzian_function_override(
    E,param
):
    return three_Lorentzian_function(E,param[0],param[1],param[2],param[3],param[4],param[5],param[6],param[7],param[8],param[9])


def two_Lorentzian_function(
    E,E1,E2,D1,D2,A1,A2,C
):
    #E is energy peak position,
    #D is the width of each peak
    #A is amplitude of each lorentzian peak
    return Lorentzian_function(E,E1,D1,A1)+Lorentzian_function(E,E2,D2,A2)+C


def Lorentzian_function(
    E,E0,D,A
):
    return A*D**2/4/((E-E0)**2+(D/2)**2)