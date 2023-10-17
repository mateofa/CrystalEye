import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad

# Define detector parameters

radius=14.4 #cm
area_det=1321 #cm2
time=100

# Define Background functions

def albedo_back(E):
    C_albedo = 1.48e-2 #Taken from Ajello 2008
    return(C_albedo/((E/33.7)**(-5)+(E/33.7)**(1.72)))

def cosmic_back(E):
    C_cosmic = 10.15e-2 #Taken from Ajello 2008
    return(C_cosmic/((E/30)**(1.32)+(E/30)**(2.88)))

#Define binning and integrate in each binning
back_int_albedo=[]
back_int_cosmic=[]

bin_center=[]
#deltaE=10.0 #Bin width in keV
#int_limits=np.linspace(10,3e4,(3e4-10)/deltaE)

int_limits_log=np.linspace(1,5,51)
int_limits = 10**(int_limits_log)
print int_limits


for i in range(np.size(int_limits)-1):
    back_int_albedo.append(time*area_det*2*np.pi*quad(albedo_back,int_limits[i],int_limits[i+1])[0])
    back_int_cosmic.append(time*area_det*2*np.pi*quad(cosmic_back,int_limits[i],int_limits[i+1])[0])

    bin_center.append(int_limits[i]+(int_limits[i+1]-int_limits[i])/2)

bin_center_MeV = np.multiply(bin_center,1e-3)

##Get total elements
back_int_albedo_total = sum(back_int_albedo)
back_int_cosmic_total = sum(back_int_cosmic)

print back_int_albedo_total
print back_int_cosmic_total


np.savetxt('Counts_Albedo_100s_1e1keV_1e5keV_50pts.txt', np.c_[bin_center_MeV,back_int_albedo],fmt='%1.10f')
np.savetxt('Counts_cosmic_100s_1e1keV_1e5keV_50pts.txt', np.c_[bin_center_MeV,back_int_cosmic],fmt='%1.8f')


plt.plot(bin_center,back_int_albedo,label='albedo',marker='o', linestyle='dashed')
plt.yscale('log')
plt.xscale('log')

#plt.savefig('Counts_Albedo_50pts.png')
