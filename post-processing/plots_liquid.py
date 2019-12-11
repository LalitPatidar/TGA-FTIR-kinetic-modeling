#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os

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
HMX = np.asarray([])
c_HONO = np.asarray([])
t_HONO = np.asarray([])
N2 = np.asarray([])
ONNO2 = np.asarray([])

with open('./mass_fractions_liquid.txt','r') as File:
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
    iHMX = header.index('HMX')
    ic_HONO = header.index('c_HONO')
    it_HONO = header.index('t_HONO')
    iN2 = header.index('N2')
    iONNO2 = header.index('ONNO2')
    
    new_lines = lines[1::2]

    for line in new_lines:
        time = np.append(time,float(line.split()[0]))
        temperature = np.append(temperature,float(line.split()[1])-273.0)
        mass_fraction = np.append(mass_fraction,100.0*float(line.split()[2]))
        mj = float(line.split()[2])*1.0
        CH2O = np.append(CH2O,mj*float(line.split()[iCH2O]))
        CO = np.append(CO,mj*float(line.split()[iCO]))
        HCN = np.append(HCN,mj*float(line.split()[iHCN]))
        N2O = np.append(N2O,mj*float(line.split()[iN2O]))
        NO = np.append(NO,mj*float(line.split()[iNO]))
        NO2 = np.append(NO2,mj*float(line.split()[iNO2]))
        H2O = np.append(H2O,mj*float(line.split()[iH2O]))
        CO2 = np.append(CO2,mj*float(line.split()[iCO2]))
        HMX = np.append(HMX,mj*float(line.split()[iHMX]))
        c_HONO = np.append(c_HONO,mj*float(line.split()[ic_HONO]))
        t_HONO = np.append(t_HONO,mj*float(line.split()[it_HONO]))
        N2 = np.append(N2,mj*float(line.split()[iN2]))
        ONNO2 = np.append(ONNO2,mj*float(line.split()[iONNO2]))

T_start = 275.0
T_end = 300.0

fig, ax1 = plt.subplots()

ax1.plot(temperature,HMX,label='HMX',marker='*',color='black')
ax1.set_ylabel('Mass fraction of HMX in liquid')
ax1.legend(loc='center left',ncol=1)
#ax1.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax1.set_xlim(T_start,T_end)
ax1.set_xlabel('Temperature($^o$C)')

ax2 = ax1.twinx()

ax2.plot(temperature,CH2O,label='CH2O',marker=">")
ax2.plot(temperature,CO,label='CO',marker="p")
ax2.plot(temperature,HCN,label='HCN',marker="<")
ax2.plot(temperature,N2O,label='N2O',marker="D")
ax2.plot(temperature,NO,label='NO',marker="^")
ax2.plot(temperature,NO2,label='NO2',marker="o")
ax2.plot(temperature,H2O,label='H2O',marker="s")
ax2.plot(temperature,CO2,label='CO2',marker="d")
ax2.plot(temperature,c_HONO,label='cis-HONO',marker="s",mfc='None')
ax2.plot(temperature,t_HONO,label='trans-HONO',marker="d",mfc='None')
ax2.plot(temperature,N2,label='N2',marker="<",mfc='None')
ax2.plot(temperature,ONNO2,label='ONNO2',marker="D",mfc='None')


#ax2.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax2.legend(loc='upper center',bbox_to_anchor=(0.5,1.2),ncol=4)
ax2.set_ylabel('Mass fraction of species in liquid')

plt.savefig('Liquids.pdf',bbox_inches='tight')

