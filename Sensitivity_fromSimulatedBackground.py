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
f_back_albedo = open('./AlbedoGamma_FluxAfterCuts.txt', 'r')
lines_back_albedo=f_back_albedo.readlines()
back_albedo_flux=[]
back_albedo_energy=[]
deltaE_albedo=[]

for x in lines_back_albedo:
    deltaE_albedo.append(float(x.split(' ')[3])-float(x.split(' ')[2]))
    #print float(x.split(' ')[4]),' ',(float(x.split(' ')[3])-float(x.split(' ')[2]))
#    back_albedo_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_albedo_flux.append(2*np.pi*float(x.split(' ')[4]))
    back_albedo_energy.append(float(x.split(' ')[1]))
#    back_albedo_energy_high.append(float(x.split(' ')[3]))
#    back_albedo_energy_low.append(float(x.split(' ')[2]))

f_back_albedo.close()

#print back_albedo_flux

## Cosmic
f_back_cosmic = open('./DiffuseGamma_FluxAfterCuts.txt', 'r')
lines_back_cosmic=f_back_cosmic.readlines()
back_cosmic_flux=[]
back_cosmic_energy=[]
deltaE_cosmic=[]

for x in lines_back_cosmic:
    deltaE_cosmic.append(float(x.split(' ')[3])-float(x.split(' ')[2]))
    #back_cosmic_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_cosmic_flux.append(2*np.pi*float(x.split(' ')[4]))

    back_cosmic_energy.append(float(x.split(' ')[1]))


f_back_cosmic.close()

## Add up each background in each energy bin
back_total_flux=[]
back_total_energy=[]
for i in range(np.size(back_cosmic_energy)) :
    back_total_flux.append(back_cosmic_flux[i]+back_albedo_flux[i])
    back_total_energy.append(back_cosmic_energy[i])


# Define Sensitivity function and create Sensitivity points in the background energy values

Sens_value=[]
Sens_value_new=[]

A_geo=1321
A_base=np.pi *  14.5**2
A_tot=A_geo + A_base
T=3.156e+7 #1 year in seconds
#T=1e+6 #Time in seconds

def Sens(e,B,deltaE):
    return  3/e * np.sqrt(B/(A_tot*T*deltaE))


#Get interpolated efficiency values in background energy bins
int_eff = np.interp(back_albedo_energy,eff_energy,eff)
int_eff_new = np.interp(back_albedo_energy,eff_energy_new,eff_new)


##Compute Sensitivity
for i in range(np.size(back_total_energy)):
    #print back_albedo_energy[i],' ',deltaE_albedo[i]
    Sens_value.append(back_albedo_energy[i]*back_albedo_energy[i]*Sens(int_eff[i],back_total_flux[i],deltaE_albedo[i]))
    Sens_value_new.append(back_albedo_energy[i]*back_albedo_energy[i]*Sens(int_eff_new[i],back_total_flux[i],deltaE_albedo[i]))


    #Sens_value.append(Sens(int_eff[i],back_int[i]))

##Get Sensitivity values for particular wavelengths
int_back_total_flux_511 = np.interp(511,back_total_energy,back_total_flux) # Background @ 511 keV
int_eff_511 = np.interp(511,eff_energy_new,eff_new) # Eff @ 511 keV
int_deltaE_511 = np.interp(511,back_total_energy,deltaE_albedo) # Background @ 511 keV

Sens_511 =  Sens(int_eff_511,int_back_total_flux_511,int_deltaE_511)
print Sens_511

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

CE_Sens_value_1MeV= np.interp(1e3,back_albedo_energy,Sens_value)*1.6e-9
AstroGam_Sens_value_1MeV= np.interp(1e3,energy_astrogam,Sens_astrogam)

ratio_CE_Astrogam = CE_Sens_value_1MeV/AstroGam_Sens_value_1MeV
#print ratio_CE_Astrogam



#Cut first point (19 to start after the drop)
Sens_value=Sens_value[19:]
Sens_value_new=Sens_value_new[19:]

back_albedo_energy=back_albedo_energy[19:]

######Sensitivity Plotting #############

fig, ax1 = plt.subplots()

#Change units
#back_energy_MeV=np.multiply(back_energy,1e-3)
#Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg=np.multiply(Sens_value,1.6e-9) #keV to erg
Sens_value_erg_new=np.multiply(Sens_value_new,1.6e-9) #keV to erg


#plt.plot(energy_astrogam,Sens_astrogam,label='eAstrogam Sensitivity 1yr ', linestyle='dashed')
#plt.plot(energy_SIP,Sens_SIP,label='SPI Sensitivity 1Ms ', linestyle='dashed')

#ax1.plot(back_albedo_energy,Sens_value_erg,label='CrystalEye Sensitivity 1yr', linestyle='dashed')
ax1.plot(back_albedo_energy,Sens_value_erg_new,label='CrystalEye Sensitivity 1yr - New Eff', linestyle='dashed')

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
