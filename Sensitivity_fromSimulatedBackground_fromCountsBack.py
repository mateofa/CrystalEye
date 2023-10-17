import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad

#Import Efficientcy Points



f = open('./Efficiency.txt', 'r')
lines=f.readlines()
eff=[]
eff_energy=[]
for x in lines:
    eff.append(float(x.split('  ')[1]))
    eff_energy.append(float(x.split('  ')[0]))
f.close()




##Import simulated and weighted background count Points - FIRST apply factor to get counts per unit cm-2 s-1

T_sim=1e2 #in seconds
A_geo=1321 # detector area


## Albedo
f_back_albedo = open('./albedo_CryEyeWeighted.txt', 'r')
lines_back_albedo=f_back_albedo.readlines()
back_albedo_count=[]
back_albedo_energy=[]
for x in lines_back_albedo:
    back_albedo_count.append(float(x.split('  ')[1])/(T_sim*A_geo))
    back_albedo_energy.append(float(x.split('  ')[0]))
f_back_albedo.close()


## Cosmic
f_back_cosmic = open('./cosmic_CryEyeWeighted.txt', 'r')
lines_back_cosmic=f_back_cosmic.readlines()
back_cosmic_count=[]
back_cosmic_energy=[]
for x in lines_back_cosmic:
    back_cosmic_count.append(1e-2*float(x.split('  ')[1])/(T_sim*A_geo))
    back_cosmic_energy.append(float(x.split('  ')[0]))
f_back_cosmic.close()

## Add up each background in each energy bin
back_total_count=[]
back_total_energy=[]
for i in range(np.size(back_cosmic_energy)) :
    back_total_count.append(back_cosmic_count[i]+back_albedo_count[i])
    back_total_energy.append(back_cosmic_energy[i])



#Define binning
int_limits_log=np.linspace(1,5,51)
int_limits = 10**(int_limits_log)
deltaE=[]
for i in range(np.size(int_limits)-1):
    deltaE.append(int_limits[i+1]-int_limits[i])

# Define Sensitivity function and create Sensitivity points in the background energy values

Sens_value=[]
A_geo=1321
T=100 #Time in seconds
def Sens(e,B,deltaE):
    return 3/(e*deltaE) * np.sqrt(B/(A_geo*T))
print deltaE

#Get interpolated efficiency values in background energy bins
#First change to keV
back_total_energy_keV = np.multiply(back_total_energy,1e3)
int_eff = np.interp(back_total_energy_keV,eff_energy,eff)

##Compute Sensitivity
for i in range(np.size(back_total_energy)):
    #print int_eff[i]
    #print back_count[i]
    #print deltaE[i]
    Sens_value.append(back_total_energy_keV[i]*back_total_energy_keV[i]*Sens(int_eff[i],back_total_count[i],deltaE[i]))
    #Sens_value.append(Sens(int_eff[i],back_int[i]))


###Plotting######

######Sensitivity Plottting #############

#Change units
#back_energy_MeV=np.multiply(back_energy,1e-3)
#Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg=np.multiply(Sens_value,1.6e-6) #MeV to erg

plt.plot(back_total_energy_keV,Sens_value_erg,label='Sensitivity 100 s', linestyle='dashed')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('E [keV]')
plt.ylabel(r' 3$\sigma$ Flux [erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
plt.legend()

plt.savefig('Sensitivity_fromWeightedBack.png')
