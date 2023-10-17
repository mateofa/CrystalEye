import numpy as np
from scipy.constants import Planck as h
from matplotlib import pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import CubicSpline


#Import Effective Area

f = open('./EffArea.dat', 'r')
lines=f.readlines()
eff=[]
eff_energy=[]
eff_energy_low=[]
eff_energy_up=[]

for x in lines:
    eff.append(float(x.split(' ')[3]))
    eff_energy.append(float(x.split(' ')[0]))
    eff_energy_low.append(float(x.split(' ')[1]))
    eff_energy_up.append(float(x.split(' ')[2]))

f.close()

#
#Import Old Sensitivity

f = open('./SensitivityCrystalEye.txt', 'r')
lines=f.readlines()
Sens_old=[]
ener_old=[]
for x in lines:
    ener_old.append(float(x.split('  ')[0]))
    Sens_old.append(float(x.split('  ')[1]))
f.close()

##Import simulated and weighted background count Points - FIRST multiply by 2pi (it's in sr-1 and we integrate over half sphere) to get counts

## Albedo
C_albedo = 1.87e-2 #Taken from Ajello 2008
def albedo_back(E):
    if(E<200):
        return(C_albedo/((E/33.7)**(-5)+(E/33.7)**(1.72)))
    if(200<=E<2e4):
        return(1.01e-4*(E/1e3)**(-1.34))
    if(2e4<=E):
         return(7.29e-4*(E/1e3)**(-2))


## Cosmic gamma
def gamma_diff_back(E):
    C_cosmic = 10.15e-2 #Taken from Ajello 2008
    return(C_cosmic/((E/30)**(1.32)+(E/30)**(2.88)))


## Cosmic gamma (LIBO)
def gamma_diff_back2(E):
	#count*cm^-2*s^1*sr^-1*keV^-1
		#cm^-2*s^1*sr^-1*keV^-1
    	if(E<=60.0):
    		f = 7.877*pow(E,-1.29)*np.exp(-E/41.13)

    	elif(E>60.0):
    		f = 4.32e-4*pow(E/60., -6.5)
    		+ 8.4e-3*pow(E/60.,-2.58)
    		+ 4.8e-4*pow(E/60.,-2.05)

    	return(f)

## Primary protons
Z=1
#e=1.602e-19
e=1
#e=1.602e-19 #(LIBO value)
F=1.23e-8
E1=5.1e5
a=0.155
#A=16.9e-7
A=2.39e-6
b=2.83
phi=6.5e5
#E2=12.25e6
E2=1.37e7 #(LIBO value)

mc2 = 938272.0813 # Proton mass in units of keV/c2
def proton_back(E):
    return(
           F*((E/1e6)**(-a))
           *(np.exp(-E/E1))**(-a+1)
           +A*((E+Z*e*phi)/1e6)**(-b)
           *((E+mc2)**2 - (mc2)**2)/((E+mc2+Z*e*phi)**2-(mc2)**2)
           *(1/(1+(E/E2)**(-12)))
           )

## Secondary protons
a_sec=0.155
E_cutSec=5.1e5
def proton_back_sec(E):
    return(
            1.23e-8*(E/1e6)**(-a_sec)*np.exp((-1)*(E/E_cutSec)**(-a_sec+1))
            )

## Primary neutrons
def neutron_back(E):
    if(10<=E<1e3):
        return(9.98e-8*(E/1e6)**(-0.5))
    if(1e3<=E<1e5):
        return(3.16e-9*(E/1e6)**(-1))
    if(1e5<=E<=1e8):
         return(3.16e-10*(E/1e6)**(-2))


##Trapped protons
f_trapped_prot=open('../../../BackgroundEnvironment/TrappedProtons.dat','r')
lines_trapped_prot=f_trapped_prot.readlines()

trapped_prot_flux=[]
trapped_prot_energy=[]

for x in lines_trapped_prot:
    trapped_prot_flux.append(float(x.split('  ')[1]))
    trapped_prot_energy.append(float(x.split('  ')[0]))
f_trapped_prot.close()

#create interpoalete funciton

natural_prot = CubicSpline(trapped_prot_energy, trapped_prot_flux, bc_type='natural')
##Trapped electrons
f_trapped_elec=open('../../../BackgroundEnvironment/TrappedElectrons.dat','r')
lines_trapped_elec=f_trapped_elec.readlines()

trapped_elec_flux=[]
trapped_elec_energy=[]

for x in lines_trapped_elec:
    trapped_elec_flux.append(float(x.split('  ')[1]))
    trapped_elec_energy.append(float(x.split('  ')[0]))
f_trapped_elec.close()

natural_elec = CubicSpline(trapped_elec_energy, trapped_elec_flux, bc_type='natural')

##Add them up
def back_total(E):
    return albedo_back(E) + gamma_diff_back2(E) + proton_back(E) + proton_back_sec(E) + neutron_back(E) + natural_elec(E) + natural_prot(E)

## Get energy bin widths from Effective Area
deltaE=[]
for i in range(np.size(eff_energy)-1):
     deltaE.append(eff_energy_up[i]-eff_energy_low[i])

# Define Sensitivity function and evaluate with proper background

Sens_value_1deg=[]
Sens_value_5deg=[]
Sens_value_90deg=[]

Sens_value_new=[]

#A_geo=1321
#A_base=np.pi *  14.5**2
#A_tot=A_geo + A_base
T=3.156e+7 #1 year in seconds
epsilon_68=0.68 #Fraction of events that are within ang resolution

theta1deg=0.017 #Angular resolution
theta5deg=0.087

deltaOmega_1deg = 2*np.pi*(1-np.cos(theta1deg))
deltaOmega_5deg = 2*np.pi*(1-np.cos(theta5deg))
deltaOmega_90deg = 2*np.pi

def Sens(B,deltaE,Aeff,deltaOmega):
    return  3/epsilon_68 * np.sqrt(B*deltaOmega/(Aeff*T*deltaE))


##Compute Sensitivity
energy_PLOT=[]
for i in range(np.size(eff_energy)):
    #print back_albedo_energy[i],' ',deltaE_albedo[i]
    if(eff[i]>0):
        energy_PLOT.append(eff_energy[i])
        Sens_value_1deg.append(eff_energy[i]*eff_energy[i]*Sens(back_total(eff_energy[i]),deltaE[i],eff[i],deltaOmega_1deg))
        Sens_value_5deg.append(eff_energy[i]*eff_energy[i]*Sens(back_total(eff_energy[i]),deltaE[i],eff[i],deltaOmega_5deg))
        Sens_value_90deg.append(eff_energy[i]*eff_energy[i]*Sens(back_total(eff_energy[i]),deltaE[i],eff[i],deltaOmega_90deg))
        #print eff[i]," ",Sens(back_total_flux[i],deltaE_albedo[i],eff[i],deltaOmega_1deg)," ",back_total_flux[i]

    #Sens_value_new.append(back_albedo_energy[i]*back_albedo_energy[i]*Sens(int_eff_new[i],back_total_flux[i],deltaE_albedo[i]))


    #Sens_value.append(Sens(int_eff[i],back_int[i]))

##Get Sensitivity values for particular wavelengths
#int_back_total_flux_511 = np.interp(511,back_total_energy,back_total_flux) # Background @ 511 keV
#int_eff_511 = np.interp(511,eff_energy_new,eff_new) # Eff @ 511 keV
#int_deltaE_511 = np.interp(511,back_total_energy,deltaE_albedo) # Background @ 511 keV

#Sens_511 =  Sens(int_eff_511,int_back_total_flux_511,int_deltaE_511)
#print Sens_511

##Add other experiments ###

## eAstrogam
f_astrogam = open('../eAstrogamData.dat', 'r')
lines_astrogam=f_astrogam.readlines()
Sens_astrogam=[]
energy_astrogam=[]
for x in lines_astrogam:
    Sens_astrogam.append(float(x.split(' ')[1]))
    energy_astrogam.append(1e3*float(x.split(' ')[0]))
f_astrogam.close()


## SIP
f_SIP = open('../SIPData.dat', 'r')
lines_SIP=f_SIP.readlines()
Sens_SIP=[]
energy_SIP=[]
for x in lines_SIP:
    Sens_SIP.append(float(x.split(' ')[1]))
    energy_SIP.append(1e3*float(x.split(' ')[0]))
f_SIP.close()

## Check ratio of values between CrysrtalEye and other experiments

#CE_Sens_value_1MeV= np.interp(1e3,back_albedo_energy,Sens_value_90deg)*1.6e-9
#AstroGam_Sens_value_1MeV= np.interp(1e3,energy_astrogam,Sens_astrogam)

#ratio_CE_Astrogam = CE_Sens_value_1MeV/AstroGam_Sens_value_1MeV
#print ratio_CE_Astrogam



#Cut first point (19 to start after the drop)
#Sens_value=Sens_value[19:]
#Sens_value_new=Sens_value_new[19:]

#back_albedo_energy=back_albedo_energy[19:]


### BACKGROUND Plotting
fig0, ax0 = plt.subplots()
    #return albedo_back(E) + gamma_diff_back2(E) + proton_back(E) + proton_back_sec(E) + neutron_back(E) + natural_elec(E) + natural_prot(E)
#ax0.plot(eff_energy,gamma_diff_back2(eff_energy))

plt.savefig('Backs.png')


######Sensitivity Plotting #############

fig, ax1 = plt.subplots()

#Change units
#back_energy_MeV=np.multiply(back_energy,1e-3)
#Sens_value_MeV=np.multiply(Sens_value,1e-3)
Sens_value_erg_1deg=np.multiply(Sens_value_1deg,1.6e-9) #keV to erg
Sens_value_erg_5deg=np.multiply(Sens_value_5deg,1.6e-9) #keV to erg
Sens_value_erg_90deg=np.multiply(Sens_value_90deg,1.6e-9) #keV to erg



Sens_value_erg_new=np.multiply(Sens_value_new,1.6e-9) #keV to erg


#plt.plot(energy_astrogam,Sens_astrogam,label='eAstrogam Sensitivity 1yr ', linestyle='dashed')
#plt.plot(energy_SIP,Sens_SIP,label='SPI Sensitivity 1Ms ', linestyle='dashed')

#ax1.plot(back_albedo_energy,Sens_value_erg,label='CrystalEye Sensitivity 1yr', linestyle='dashed')
ax1.plot(energy_PLOT,Sens_value_erg_1deg,label=r"CE Sensitivity 1yr - $\Delta \theta$= 1 deg", linestyle='dashed')
ax1.plot(energy_PLOT,Sens_value_erg_5deg,label=r"CE Sensitivity 1yr - $\Delta \theta$= 5 deg", linestyle='dashed')
ax1.plot(energy_PLOT,Sens_value_erg_90deg,label=r"CE Sensitivity 1yr - $\Delta \theta$= 90 deg (Full Sky)", linestyle='dashed')
ax1.plot(ener_old,Sens_old,label=r"Previous CE Sensitivity 1yr - (Full Sky)", linestyle='solid')

ax1.set_title('$\epsilon_{68}=0.68$')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('E [keV]')
ax1.set_ylabel(r' 3$\sigma$ Sensitivity [erg cm$^{-2}$ s$^{-1}$]')
ax1.legend(loc=2, prop={'size': 10})
#ax1.legend()
plt.xlim([5,1e5 ])
plt.ylim([1e-14,1e-9 ])


#ax2 = ax1.twinx()
#ax2.plot(eff_energy,eff,label='Efficiency', linestyle='dashed',color='red')
#ax2.set_yscale('log')
#ax2.set_xscale('log')
#ax2.set_xlabel('E [keV]')
#ax2.set_ylabel('Efficiency')
#ax2.legend()



plt.savefig('Sensitivity_fromSimulatedBack.png')
