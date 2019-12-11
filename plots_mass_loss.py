#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os

time = np.asarray([])
#temperature = np.asarray([])
#mass_fraction = np.asarray([])

temperature = np.linspace(250.0,275.0,20)
mass_fraction = np.linspace(100.,100.,20)

CH2O = np.asarray([])
CO = np.asarray([])
HCN = np.asarray([])
N2O = np.asarray([])
NO = np.asarray([])
NO2 = np.asarray([])
H2O = np.asarray([])
CO2 = np.asarray([])

with open('./output_TGA/mole_fractions_gas.txt','r') as File:
    lines = File.readlines()
    
    header = lines[0].split()
    
    iCH2O = header.index('CH2O')
    iCO = header.index('CO')
    iHCN = header.index('HCN')
    iN2O = header.index('N2O')
    iNO = header.index('NO')
    iNO2 = header.index('NO2')
    iH2O = header.index('H2O')
    iCO2 = header.index('CO2')
    
    new_lines = lines[1::3]
    #temp = lines[70::100]
    
    #new_lines.extend(temp)

    for line in new_lines:
        time = np.append(time,float(line.split()[0]))
        temperature = np.append(temperature,float(line.split()[1])-273.0)
        mass_fraction = np.append(mass_fraction,100.0*float(line.split()[2]))
        CH2O = np.append(CH2O,float(line.split()[iCH2O]))
        CO = np.append(CO,float(line.split()[iCO]))
        HCN = np.append(HCN,float(line.split()[iHCN]))
        N2O = np.append(N2O,float(line.split()[iN2O]))
        NO = np.append(NO,float(line.split()[iNO]))
        NO2 = np.append(NO2,float(line.split()[iNO2]))
        H2O = np.append(H2O,float(line.split()[iH2O]))
        CO2 = np.append(CO2,float(line.split()[iCO2]))

HMXc = np.asarray([])
Tc = np.asarray([])
with open('./output_TGA/mass_fractions_liquid.txt','r') as File:
    lines = File.readlines()
    
    header = lines[0].split()
    
    iHMX = header.index('HMX')
    for line in lines[1:]:
        Tc = np.append(Tc,float(line.split()[1])-273.0)
        HMXc = np.append(HMXc,float(line.split()[iHMX])*100.0)

#---------------------------------------------------------------------------------------------
Texp = np.asarray([])
m1 = np.asarray([])
m2 = np.asarray([])
m3 = np.asarray([])
m4 = np.asarray([])
m5 = np.asarray([])
m6 = np.asarray([])
m7 = np.asarray([])
m8 = np.asarray([])
m9 = np.asarray([])

dsc1 = np.asarray([])
dsc2 = np.asarray([])
dsc3 = np.asarray([])
dsc4 = np.asarray([])
dsc5 = np.asarray([])
dsc6 = np.asarray([])
dsc7 = np.asarray([])
dsc8 = np.asarray([])
dsc9 = np.asarray([])

start = 1
step = 1

path_mass_loss = './mass_loss_data'

with open(os.path.join(path_mass_loss,'HMX_5KPM_1198ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        Texp = np.append(Texp,float(line.split()[0]))
        dsc1 = np.append(dsc1,float(line.split()[2]))
        m1 = np.append(m1,float(line.split()[3])-15.0)

with open(os.path.join(path_mass_loss,'HMX_5KPM_1180ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc2 = np.append(dsc2,float(line.split()[2]))
        m2 = np.append(m2,float(line.split()[3])-1.0)
        
with open(os.path.join(path_mass_loss,'HMX_5KPM_1104ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc3 = np.append(dsc3,float(line.split()[2]))
        m3 = np.append(m3,float(line.split()[3])-8.0)

with open(os.path.join(path_mass_loss,'HMX_10KPM_1050ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc4 = np.append(dsc4,float(line.split()[2])-2.5)
        m4 = np.append(m4,float(line.split()[3]))

with open(os.path.join(path_mass_loss,'HMX_10KPM_1024ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc5 = np.append(dsc5,float(line.split()[2])-2.2)
        m5 = np.append(m5,float(line.split()[3])-7.0)
        
with open(os.path.join(path_mass_loss,'HMX_10KPM_1094ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc6 = np.append(dsc6,float(line.split()[2])-2.1)
        m6 = np.append(m6,float(line.split()[3])-6.0)
        
#with open(os.path.join(path_mass_loss,'HMX_15KPM_1176ug.txt'),'r') as File:
#    lines = File.readlines()
#    lines = lines[start::step]
#    for line in lines:
#        dsc7 = np.append(dsc7,float(line.split()[2]))
#        m7 = np.append(m7,float(line.split()[3])+3.0)

with open(os.path.join(path_mass_loss,'HMX_15KPM_1158ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc8 = np.append(dsc8,float(line.split()[2])-7.2)
        m8 = np.append(m8,float(line.split()[3])+1.0)
        
with open(os.path.join(path_mass_loss,'HMX_15KPM_1142ug.txt'),'r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        dsc9 = np.append(dsc9,float(line.split()[2])-7.2)
        m9 = np.append(m9,float(line.split()[3])+3.0)

m5_average = np.mean([m1,m2,m3],axis=0)
m10_average = np.mean([m4,m5,m6],axis=0)
m15_average = np.mean([m8,m9],axis=0)

dsc5_average = np.mean([dsc1,dsc2,dsc3],axis=0)
dsc10_average = np.mean([dsc4,dsc5,dsc6],axis=0)
dsc15_average = np.mean([dsc8,dsc9],axis=0)
#---------------------------------------------------------------------------------------------

T_start = 250.0
T_end = 300.0

fig, ax1 = plt.subplots()
"""
ax1.plot(temperature,100.0*CH2O,label='CH2O',marker=">")
ax1.plot(temperature,100.0*CO,label='CO',marker="p")
ax1.plot(temperature,100.0*HCN,label='HCN',marker="<")
ax1.plot(temperature,100.0*N2O,label='N2O',marker="D")
ax1.plot(temperature,100.0*NO,label='NO',marker="^")
ax1.plot(temperature,100.0*NO2,label='NO2',marker="o")
ax1.plot(temperature,100.0*H2O,label='H2O',marker="s")
ax1.plot(temperature,100.0*CO2,label='CO2',marker="d")
"""
ax1.plot(Texp[::4],m15_average[::4],label='Experimental',marker='*',color='black',linestyle='None')
ax1.plot(temperature,mass_fraction,label='Computational',color='black')

#ax1.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax1.set_xlim(T_start,T_end)
ax1.set_xlabel('Temperature($^o$C)')
ax1.set_ylabel('Condensed-phase mass (%)')
ax1.legend(loc='upper center',bbox_to_anchor=(0.5,1.15),ncol=4)

#ax2 = ax1.twinx()
#ax2.plot(temperature,mass_fraction,label='TGA',color='black')

#ax2.plot(Texp,m5_average,label='TGA',marker="*",color='black')
#ax2.plot(Texp,m10_average,label='TGA',marker="*",color='black')
#ax2.plot(Texp[::4],m15_average[::4],label='TGA',marker='*',color='black')

#ax2.plot(Tc,HMXc,color='magenta')
#ax2.legend()
#ax2.set_ylabel('Condensed-phase mass (%)')
plt.savefig('./output_TGA/TG.pdf')

