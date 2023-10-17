import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad

#Import Efficientcy Points

f = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/EfficiencyFiles/EfficiencyAll_ContinuousSpectrum/Efficiency_0_0.txt', 'r')
lines=f.readlines()
eff=[]
eff_energy=[]
for x in lines:
    eff.append(float(x.split('  ')[1]))
    eff_energy.append(float(x.split('  ')[0]))
f.close()

#Import eAstrogam background points

f_astroGam = open('Back_eAstrogam.dat', 'r')
lines_astroGam=f_astroGam.readlines()
energy_astroGam=[]
back_astroGam=[]
for x in lines_astroGam:
    energy_astroGam.append(1e3*float(x.split(' ')[0]))
    back_astroGam.append(1e-7*float(x.split(' ')[1])/(2*np.pi))
f_astroGam.close()

# Define Background functions

def albedo_back(E):
    C_albedo = 1.48e-2 #Taken from Ajello 2008
    return(C_albedo/((E/33.7)**(-5)+(E/33.7)**(1.72)))

def cosmic_back(E):
    C_cosmic = 10.15e-2 #Taken from Ajello 2008
    return(C_cosmic/((E/30)**(1.32)+(E/30)**(2.88)))


def back_total(E):
    return cosmic_back(E) + albedo_back(E)


#Define binning and integrate in each binning
back_int=[]
bin_center=[]

deltaE=10.0 #Bin width in keV
int_limits=np.linspace(2,3e4,(3e4-10)/deltaE)


for i in range(np.size(int_limits)-1):
    back_int.append(2*np.pi*quad(back_total,int_limits[i],int_limits[i+1])[0])
    bin_center.append(int_limits[i]+(int_limits[i+1]-int_limits[i])/2)


#Get interpolated efficiency values in bin_centers
int_eff = np.interp(bin_center,eff_energy,eff)

# Define Sensitivity function and create Sensitivity points in the bin_centers
Sens_value=[]
A_geo=1321
T=1e+6 #Time in seconds
def Sens(e,B):
    return 3/(e*deltaE) * np.sqrt(B/(A_geo*T))

for i in range(np.size(bin_center)):
    Sens_value.append(bin_center[i]*bin_center[i]*Sens(int_eff[i],back_int[i]))
    #Sens_value.append(Sens(int_eff[i],back_int[i]))


###Plotting######

######Sensitivity Plottting #############

#Change units
bin_center_MeV=np.multiply(bin_center,1e-3)
Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg=np.multiply(Sens_value,1.6e-9)

plt.plot(bin_center_MeV,Sens_value_erg,label='Sensitivity 1e+6 s', linestyle='dashed')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('E [MeV]')
plt.ylabel(r' 3$\sigma$ Flux [erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
plt.legend()

plt.savefig('Sensitivity_Continuous.png')

###Backgound plotting###########

#Define energy points and get log-flux

E_full = np.logspace(1,4,100)
albedo_flux = albedo_back(E_full)
cosmic_flux =cosmic_back(E_full)
total_flux =back_total(E_full)

#log_E = np.log10(E_full)
#log_albedo_flux = np.log10(albedo_back(E_full))
#log_cosmic_flux = np.log10(cosmic_back(E_full))

plt.plot(energy_astroGam,back_astroGam,label='Cosmic eAstrogam', linestyle='dashed')

plt.plot(E_full,albedo_flux,label='albedo', linestyle='dashed')
plt.plot(E_full,cosmic_flux,label='cosmic', linestyle='dashed')
plt.plot(E_full,total_flux,label='total')
#plt.plot(bin_center,integral,label='integral')
plt.legend()
plt.xlabel('E [keV]')
plt.ylabel(r'Flux [ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$ sr$^{-1}$]')

plt.yscale('log')
plt.xscale('log')

plt.savefig('Flux_Back.png')
