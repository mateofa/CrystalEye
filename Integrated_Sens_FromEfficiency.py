import glob,os
import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad

angle=0

f = open('./EfficiencyTXT2/Efficiency_Angle_0.txt', 'r')
lines=f.readlines()
eff=[]
eff_energy=[]
for x in lines:
    eff.append(float(x.split('  ')[1]))
    eff_energy.append(float(x.split('  ')[0]))
f.close()




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



A_geo=1321
A_tot=A_geo
#T=3.156e+7 #1 year in seconds
T=1 #Time in seconds

def Sens(e,B,deltaE):
    return  3/e * np.sqrt(B/(A_tot*T*deltaE))


#Get interpolated efficiency values in background energy bins
int_eff = np.interp(back_albedo_energy,eff_energy,eff)

##Compute Sensitivity
Sens_value_int_20keV=0
Sens_value_int_Tot=0

Sens_value_20keV=[]
E_int_20keV=[]

Sens_value_Tot=[]
E_int_Tot=[]

Emin=0
Emax=20
for i in range(np.size(back_total_energy)):
    #print back_albedo_energy[i],' ',deltaE_albedo[i]
    Sens_value.append(back_albedo_energy[i]*back_albedo_energy[i]*Sens(int_eff[i],back_total_flux[i],deltaE_albedo[i]))
    if (Emin<back_albedo_energy[i]<Emax and (int_eff[i] != 0)):
        E_int_20keV.append(back_albedo_energy[i])
        Sens_value_20keV.append(Sens(int_eff[i],back_total_flux[i],deltaE_albedo[i]))
        Sens_value_int_20keV += Sens(int_eff[i],back_total_flux[i],deltaE_albedo[i])*deltaE_albedo[i]

    if ((int_eff[i] != 0)):
        E_int_Tot.append(back_albedo_energy[i])
        Sens_value_Tot.append(Sens(int_eff[i],back_total_flux[i],deltaE_albedo[i]))
        Sens_value_int_Tot += Sens(int_eff[i],back_total_flux[i],deltaE_albedo[i])*deltaE_albedo[i]
        #Integrate in FWHF energies_hawc


print Sens_value_int_20keV

##Get integral using trapezoid Method
Int_20=np.trapz(Sens_value_20keV,x=E_int_20keV)
print Int_20


#print angle," ",Sens_value_int_20keV," ",Sens_value_int_Tot
#f_results=open("Integrated_Sensitivity_round2.txt", 'a+')
#print >> f_results, angle," ",Sens_value_int_20keV," ",Sens_value_int_Tot
#f_results.close()

f_sens=open("Integrated_Sensitivity.txt", 'r')
lines_sens=f_sens.readlines()[1:]
angle=[]
sens_20=[]
sens_tot=[]
for x in lines_sens:
    angle.append(float(x.split('  ')[0]))
    sens_20.append(float(x.split('  ')[1]))
    sens_tot.append(float(x.split('  ')[2]))

f_sens.close()
plt.title("Integrated Sens. 1sec Exposure")
plt.plot(angle,sens_20,label='Integrated Sensitivity 10 keV - 20 keV', linestyle='dashed')
plt.plot(angle,sens_tot,label='Integrated Sensitivity 10 keV - 100 MeV', linestyle='dashed')
plt.xlabel('angle [deg]')
plt.ylabel(r' Int. Sens. [ph cm$^{-2}$ s$^{-1}$]')
plt.legend()

plt.savefig('SensVSanlge.png')
