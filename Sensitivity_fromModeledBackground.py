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

#f_new = open('./EfficiencyTXT2/Efficiency_Angle_0.txt', 'r')
f_new = open('/Users/mateo/Documents/GSSI/Nuses/CrystalEye/EfficiencyFiles/EfficiencyAll_ContinuousSpectrum/Efficiency_0_0.txt', 'r')

lines_new=f_new.readlines()
eff_new=[]
eff_energy_new=[]
for x in lines_new:
    eff_new.append(float(x.split('  ')[1]))
    eff_energy_new.append(float(x.split('  ')[0]))
f_new.close()



##Import simulated and weighted background count Points - FIRST multiply by 2pi (it's in sr-1 and we integrate over half sphere) to get counts

## Albedo
def albedo_back(E):
    C_albedo = 1.48e-2 #Taken from Ajello 2008
    return(C_albedo/((E/33.7)**(-5)+(E/33.7)**(1.72)))
#print back_albedo_flux

## Cosmic
def cosmic_back(E):
    C_cosmic = 10.15e-2 #Taken from Ajello 2008
    return(C_cosmic/((E/30)**(1.32)+(E/30)**(2.88)))

##Add them up
def back_total(E):
    return cosmic_back(E) + albedo_back(E)



# Define Sensitivity function and create Sensitivity points in the background energy values

Sens_value=[]
Sens_value_new=[]

A_geo=1321
A_base=np.pi *  14.5**2
A_tot=A_geo + A_base
T=3.156e+7 #1 year in seconds
#T=1e+6 #Time in seconds


##Define binning
bin_limits=np.logspace(1,4,51)
bin_center_array=[]
bin_width_array=[]

for i in range(len(bin_limits)-1):
    bin_center=bin_limits[i+1]-bin_limits[i]/2
    bin_width=bin_limits[i+1]-bin_limits[i]
    bin_width_array.append(bin_width)
    bin_center_array.append(bin_center)


def Sens(e,B,deltaE):
    return  3/e * np.sqrt(B/(A_tot*T*deltaE))


#Get interpolated efficiency values in background energy bins
int_eff = np.interp(bin_center_array,eff_energy,eff)
int_eff_new = np.interp(bin_center_array,eff_energy_new,eff_new)



##Compute Sensitivity
for i in range(len(bin_center_array)):
    #print back_albedo_energy[i],' ',deltaE_albedo[i]
    Sens_value.append(bin_center_array[i]*bin_center_array[i]*Sens(int_eff_new[i],back_total(bin_center_array[i]),bin_width_array[i]))


    #Sens_value.append(Sens(int_eff[i],back_int[i]))

##Get Sensitivity values for particular wavelengths
int_eff_511 = np.interp(511,eff_energy_new,eff_new) # Eff @ 511 keV
int_deltaE_511 = np.interp(511,bin_center_array,bin_width_array) # Background @ 511 keV

Sens_511 =  Sens(int_eff_511,back_total(511),int_deltaE_511)

##Add other experiments ###

## eAstrogam
f_astrogam = open('./eAstrogamData.dat', 'r')
lines_astrogam=f_astrogam.readlines()
Sens_astrogam=[]
energy_astrogam=[]
for x in lines_astrogam:
    Sens_astrogam.append(float(x.split(' ')[1]))
    energy_astrogam.append(1e3*float(x.split(' ')[0]))
f_astrogam.close()


## SIP
f_SIP = open('./SIPData.dat', 'r')
lines_SIP=f_SIP.readlines()
Sens_SIP=[]
energy_SIP=[]
for x in lines_SIP:
    Sens_SIP.append(float(x.split(' ')[1]))
    energy_SIP.append(1e3*float(x.split(' ')[0]))
f_SIP.close()

## Check ratio of values between CrysrtalEye and other experiments

CE_Sens_value_1MeV= np.interp(1e3,bin_center_array,Sens_value)*1.6e-9
AstroGam_Sens_value_1MeV= np.interp(1e3,energy_astrogam,Sens_astrogam)

ratio_CE_Astrogam = CE_Sens_value_1MeV/AstroGam_Sens_value_1MeV
print ratio_CE_Astrogam



#Cut first point (19 to start after the drop)
Sens_value=Sens_value[19:]

bin_center_array=bin_center_array[19:]

######Sensitivity Plotting #############

fig, ax1 = plt.subplots()

#Change units
#back_energy_MeV=np.multiply(back_energy,1e-3)
#Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg=np.multiply(Sens_value,1.6e-9) #keV to erg


#plt.plot(energy_astrogam,Sens_astrogam,label='eAstrogam Sensitivity 1yr ', linestyle='dashed')
#plt.plot(energy_SIP,Sens_SIP,label='SPI Sensitivity 1Ms ', linestyle='dashed')

#ax1.plot(back_albedo_energy,Sens_value_erg,label='CrystalEye Sensitivity 1yr', linestyle='dashed')
ax1.plot(bin_center_array,Sens_value_erg,label='CrystalEye Sensitivity 1yr - New Eff', linestyle='dashed')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('E [keV]')
ax1.set_ylabel(r' 3$\sigma$ Sensitivity [erg cm$^{-2}$ s$^{-1}$]')
ax1.legend()
plt.xlim([10,1e5 ])

#ax2 = ax1.twinx()
#ax2.plot(eff_energy,eff,label='Efficiency', linestyle='dashed',color='red')
#ax2.set_yscale('log')
#ax2.set_xscale('log')
#ax2.set_xlabel('E [keV]')
#ax2.set_ylabel('Efficiency')
#ax2.legend()



plt.savefig('Sensitivity_fromWeightedBack.png')
