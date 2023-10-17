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
T_year=3.156e+7 #1 year in seconds
T_mega=1e+6 #Time in seconds


def Sens(e,B,deltaE):
    return  3/e * np.sqrt(B/(A_tot*T_mega*deltaE))


#Get interpolated values for efficiency, background and binning in new range
new_range=np.logspace(-1,4.65,1e6)

##get bin width

bin_width_new=[(new_range[0]-0.1)/2 + (new_range[1]-new_range[0])/2] ## First bin width
for i in range(np.size(new_range)-2):
    bin_width_new.append((new_range[i+1]-new_range[i])/2+(new_range[i+2]-new_range[i+1])/2)
bin_width_new.append((10**4.65-new_range[-1])/2 + (new_range[-1]-new_range[-2])/2) ## last bin width


##get bin width Second Method
bin_width_new_b=[x - new_range[i - 1] for i, x in enumerate(new_range)][1:]
#print bin_width_new[1],' ',bin_width_new_b[1]

#print bin_width_new
int_eff = np.interp(new_range,eff_energy,eff)
int_back_total_flux = np.interp(new_range,back_total_energy,back_total_flux)
int_deltaE = np.interp(new_range,back_total_energy,deltaE_albedo)

##Compute Sensitivity
Sens_value_int_20keV=0
Sens_value_int_Tot=0
Sens_line=0

##Define iintegral range for line Sensitivity
line = 4438
FWHM = line*0.13
E_low = line - FWHM
E_high =line + FWHM

Sens_value=[]
Sens_value_line=[]
energy_line=[]

for i in range(np.size(new_range)):

    Sens_value.append(Sens(int_eff[i],int_back_total_flux[i],int_deltaE[i]))
    if(E_low<=new_range[i]<=E_high):
        #print Sens(int_eff[i],int_back_total_flux[i],bin_width_new[i]),' ',bin_width_new[i]
        energy_line.append(new_range[i])
        Sens_value_line.append(Sens(int_eff[i],int_back_total_flux[i],int_deltaE[i]))
        #Sens_line += Sens(int_eff[i],int_back_total_flux[i],int_deltaE[i])*bin_width_new[i]
        Sens_line += Sens(int_eff[i],int_back_total_flux[i],bin_width_new_b[i])*bin_width_new_b[i]
        #Sens_line += Sens(int_eff[i],int_back_total_flux[i],1)





##Try numpy trapz for integral
Sens_line_2=np.trapz(Sens_value_line,x=energy_line)
print Sens_line_2
MeVtoERG=1.60218e-6
ERGtoKeV=6.242e+8
##Compare with eAstrogam Sensitivity


## eAstrogam
f_astrogam = open('./eAstrogamData.dat', 'r')
lines_astrogam=f_astrogam.readlines()
Sens_astrogam=[]
energy_astrogam=[]
for x in lines_astrogam:
        Sens_astrogam.append(float(x.split(' ')[1])/(new_range[i]*new_range[i]))
        energy_astrogam.append(1e3*float(x.split(' ')[0]))
f_astrogam.close()

#Interpolate eAstrogam
interp_Astrogam = np.interp(new_range,energy_astrogam,Sens_astrogam)

T_factor=np.sqrt(T_year/T_mega)

Sens_astrogam_line=[]
energy_astrogam_line=[]
for i in range(np.size(new_range)):
    if(E_low<new_range[i]<E_high):
        Sens_astrogam_line.append(T_factor*ERGtoKeV*float(x.split(' ')[1])/(new_range[i]*new_range[i]))
        energy_astrogam_line.append(new_range[i])
f_astrogam.close()

##Check Integral for eAstrogam
Sens_line_eAstro=np.trapz(Sens_astrogam_line,x=energy_astrogam_line)
print Sens_line_eAstro

print Sens_line_2/Sens_line_eAstro




#print angle," ",Sens_value_int_20keV," ",Sens_value_int_Tot
#f_results=open("Integrated_Line.txt", 'a+')
#print >> f_results, angle," ",Sens_value_int_20keV," ",Sens_value_int_Tot
#f_results.close()

##Plot Sensitivity
fig, ax1 = plt.subplots()
ax1.set_yscale('log')
ax1.set_xscale('log')

#plt.scatter(new_range,Sens_value,label='Sensitivity')
plt.scatter(energy_line,Sens_value_line,label='Sensitivity CE' )
plt.scatter(energy_astrogam_line,Sens_astrogam_line,label='Sens Astrogram', marker='o')
#plt.plot(back_total_energy,deltaE_albedo,label='binning', marker='o')

plt.savefig('Sens_Test.png')
