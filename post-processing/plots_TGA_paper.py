#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import AutoMinorLocator


time = np.asarray([])
temperature = np.asarray([])
mass_fraction = np.asarray([])
CH2O = np.asarray([])
CO = np.asarray([])
HCN = np.asarray([])
N2O = np.asarray([])
NO = np.asarray([])
NO2 = np.asarray([])
H2O = np.asarray([])
CO2 = np.asarray([])

with open('./mole_fractions_gas.txt','r') as File:
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
    
    new_lines = lines[1::2]
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

T_start = 275.0
T_end = 300.0

fig, ax1 = plt.subplots()
ax1.xaxis.set_minor_locator(AutoMinorLocator())

ax1.plot(temperature,100.0*CH2O,label='CH$_2$O',marker=">")
ax1.plot(temperature,100.0*CO,label='CO',marker="p")
ax1.plot(temperature,100.0*HCN,label='HCN',marker="<")
ax1.plot(temperature,100.0*N2O,label='N$_2$O',marker="D")
ax1.plot(temperature,100.0*NO,label='NO',marker="^")
ax1.plot(temperature,100.0*NO2,label='NO$_2$',marker="o")
ax1.plot(temperature,100.0*H2O,label='H$_2$O',marker="s")
ax1.plot(temperature,100.0*CO2,label='CO$_2$',marker="d")

ax1.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax1.set_xlim(T_start,T_end)
ax1.set_xlabel('Temperature($^o$C)')
ax1.set_ylabel('Mole fraction (%)')
ax1.legend(loc='upper center',bbox_to_anchor=(0.5,1.15),ncol=4)

ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(direction='in',which='both')

textstr = '\n'.join(['','Computational'])
ax1.text(0.75,0.92,textstr,transform=ax1.transAxes)

plt.savefig('Gases.pdf')

