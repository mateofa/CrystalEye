    import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad


##Import generated background count from model

## albedo_model
f_back_albedo_model = open('./Counts_albedo_100s_1e1keV_1e5keV_50pts.txt', 'r')
lines_back_albedo_model=f_back_albedo_model.readlines()
back_albedo_model_count=[]
back_albedo_model_energy=[]
for x in lines_back_albedo_model:
    back_albedo_model_count.append(float(x.split(' ')[1]))
    back_albedo_model_energy.append(float(x.split(' ')[0])*1e3) #Import energy in keV
f_back_albedo_model.close()


## cosmic_model
f_back_cosmic_model = open('./Counts_cosmic_100s_1e1keV_1e5keV_50pts.txt', 'r')
lines_back_cosmic_model=f_back_cosmic_model.readlines()
back_cosmic_model_count=[]
back_cosmic_model_energy=[]
for x in lines_back_cosmic_model:
    back_cosmic_model_count.append(float(x.split(' ')[1]))
    back_cosmic_model_energy.append(float(x.split(' ')[0])*1e3) #Import energy in keV
f_back_cosmic_model.close()

## Add up each background in each energy bin
back_total_model_count=[]
back_total_model_energy=[]
for i in range(np.size(back_cosmic_model_energy)) :
    back_total_model_count.append(back_cosmic_model_count[i]+back_albedo_model_count[i])
    back_total_model_energy.append(back_cosmic_model_energy[i])


##Import simulated and weighted background count Points

## albedo
f_back_albedo = open('./albedo_CryEyeWeighted.txt', 'r')
lines_back_albedo=f_back_albedo.readlines()
back_albedo_count=[]
back_albedo_energy=[]
for x in lines_back_albedo:
    back_albedo_count.append(float(x.split('  ')[1]))
    back_albedo_energy.append(float(x.split('  ')[0])*1e3) #Import energy in keV
f_back_albedo.close()


## Cosmic
f_back_cosmic = open('./cosmic_CryEyeWeighted.txt', 'r')
lines_back_cosmic=f_back_cosmic.readlines()
back_cosmic_count=[]
back_cosmic_energy=[]
for x in lines_back_cosmic:
    back_cosmic_count.append(float(x.split('  ')[1]))
    back_cosmic_energy.append(float(x.split('  ')[0])*1e3) #Import energy in keV
f_back_cosmic.close()

## Add up each background in each energy bin
back_total_count=[]
back_total_energy=[]
for i in range(np.size(back_cosmic_energy)) :
    back_total_count.append(back_cosmic_count[i]+back_albedo_model_count[i])
    back_total_energy.append(back_cosmic_energy[i])


###Plotting######

###Backgound plotting###########

#Define energy points and get log-flux

plt.plot(back_albedo_model_energy,back_albedo_model_count,label='albedo from Model', linestyle='solid')
plt.plot(back_cosmic_model_energy,back_cosmic_model_count,label='cosmic from Model', linestyle='solid')

plt.plot(back_albedo_energy,back_albedo_count,label='albedo after cuts', linestyle='dashed')
plt.plot(back_cosmic_energy,back_cosmic_count,label='cosmic after cuts', linestyle='dashed')

plt.yscale('log')
plt.xscale('log')


#plt.plot(bin_center,integral,label='integral')
plt.legend()
plt.xlabel('E [MeV]')
plt.ylabel('Counts (100s)')

#plt.ylabel(r'Flux [ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$ sr$^{-1}$]')

plt.yscale('log')
plt.xscale('log')

plt.savefig('BackgroundCounts_Compared.png')
