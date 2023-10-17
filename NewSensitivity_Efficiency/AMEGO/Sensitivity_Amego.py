from scipy.integrate import quad
from sympy import symbols, solve
import numpy as np


### This script calculates Sensitivity as it is explained in Martinez-Castellanos ApJ 2022, AMEGO PAPER

##COMPTOM=NIZE MODEL FOR GRB prefered over Band model

A = 1.48e-2 # s-1 cm-2 between 10-1000 keV
alpha = -0.37
Epeak = 636

def GRB(x):
     return A*(x/100)**(alpha)*np.exp(-(alpha+2)*x/Epeak)


# First integrate between 100 keV and 1 MeV and over time
Tobs=1
F_GRB = quad(GRB, 100, 1000) * Tobs
F_GRB=F_GRB[0]


##EFFECTIVE AREA

f = open('../EffArea.dat', 'r')
lines=f.readlines()
effA=[]
effA_energy=[]

Emin_sub = 100
Emax_sub = 1000

effA_sub=[]
effA_energy_sub=[]

for x in lines:
    effA.append(float(x.split(' ')[3]))
    effA_energy.append(float(x.split(' ')[0]))
    if(Emin_sub<float(x.split(' ')[0])<Emax_sub):
        effA_sub.append(float(x.split(' ')[3]))
        effA_energy_sub.append(float(x.split(' ')[0]))
f.close()

#Spectrum averaged effective area
effA_avg=np.trapz(effA_sub,x=effA_energy_sub)/(Emax_sub-Emin_sub)
print effA_avg



##IMPORT BACKGROUNDS
A_factor_upper=(2*np.pi*14.5**2)*(2*np.pi)
A_factor_lower=(np.pi*14.5**2)*(2*np.pi)

## Albedo
f_back_albedo = open('../bkg_reweight/AlbedoGamma.dat', 'r')
lines_back_albedo=f_back_albedo.readlines()
back_albedo_flux=[]
back_albedo_energy=[]
deltaE_albedo=[]

for x in lines_back_albedo:
    deltaE_albedo.append(float(x.split(' ')[2])-float(x.split(' ')[1]))
    #print float(x.split(' ')[4]),' ',(float(x.split(' ')[3])-float(x.split(' ')[2]))
#    back_albedo_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_albedo_flux.append(float(x.split(' ')[3])*A_factor_lower)
    back_albedo_energy.append(float(x.split(' ')[0]))
#    back_albedo_energy_high.append(float(x.split(' ')[3]))
#    back_albedo_energy_low.append(float(x.split(' ')[2]))

f_back_albedo.close()

#print back_albedo_flux

## Cosmic gamma
f_back_cosmic = open('../bkg_reweight/DiffuseGamma.dat', 'r')
lines_back_cosmic=f_back_cosmic.readlines()
back_cosmic_flux=[]
back_cosmic_energy=[]
deltaE_cosmic=[]

for x in lines_back_cosmic:
    deltaE_cosmic.append(float(x.split(' ')[2])-float(x.split(' ')[1]))
    #back_cosmic_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_cosmic_flux.append(float(x.split(' ')[3])*A_factor_upper)
    back_cosmic_energy.append(float(x.split(' ')[0]))
f_back_cosmic.close()

## Primary protons
f_back_protons = open('../bkg_reweight/PrimaryProton.dat', 'r')
lines_back_protons=f_back_protons.readlines()
back_protons_flux=[]
back_protons_energy=[]
deltaE_protons=[]

for x in lines_back_protons:
    deltaE_protons.append(float(x.split(' ')[2])-float(x.split(' ')[1]))
    #back_protons_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_protons_flux.append(float(x.split(' ')[3])*A_factor_upper)
    back_protons_energy.append(float(x.split(' ')[0]))
f_back_protons.close()

## Primary protons_sec
f_back_protons_sec = open('../bkg_reweight/SecondaryProton.dat', 'r')
lines_back_protons_sec=f_back_protons_sec.readlines()
back_protons_sec_flux=[]
back_protons_sec_energy=[]
deltaE_protons_sec=[]

for x in lines_back_protons_sec:
    deltaE_protons_sec.append(float(x.split(' ')[2])-float(x.split(' ')[1]))
    #back_protons_sec_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_protons_sec_flux.append(float(x.split(' ')[3])*A_factor_lower)
    back_protons_sec_energy.append(float(x.split(' ')[0]))
f_back_protons_sec.close()

## Primary neutrons
f_back_neutrons = open('../bkg_reweight/Neutron.dat', 'r')
lines_back_neutrons=f_back_neutrons.readlines()
back_neutrons_flux=[]
back_neutrons_energy=[]
deltaE_neutrons=[]

for x in lines_back_neutrons:
    deltaE_neutrons.append(float(x.split(' ')[2])-float(x.split(' ')[1]))
    #back_neutrons_flux.append(2*np.pi*float(x.split(' ')[4])*(float(x.split(' ')[3])-float(x.split(' ')[2])))
    back_neutrons_flux.append(float(x.split(' ')[3])*A_factor_lower)
    back_neutrons_energy.append(float(x.split(' ')[0]))
f_back_neutrons.close()


## Add up all backgrounds in each energy bin, note that these come in ph keV-1 s-1

back_total_flux=[]
back_total_energy=[]

back_total_energy_sub=[]
back_total_flux_sub=[]

for i in range(np.size(back_cosmic_energy)) :
    back_total_flux.append(back_cosmic_flux[i]+back_albedo_flux[i]+back_protons_flux[i]+back_protons_sec_flux[i]+back_neutrons_flux[i])
    back_total_energy.append(back_cosmic_energy[i])
    if(Emin_sub<back_cosmic_energy[i]<Emax_sub):
        back_total_flux_sub.append(back_cosmic_flux[i]+back_albedo_flux[i]+back_protons_flux[i]+back_protons_sec_flux[i]+back_neutrons_flux[i])
        back_total_energy_sub.append(back_cosmic_energy[i])

## Integrate over energy range to get the RATE

back_rate =np.trapz(back_total_flux_sub,x=back_total_energy_sub)
print  back_total_Int



##NOW LET'S SET THE EQUATION

S=F*effA_avg*T
B=back_rate*T
