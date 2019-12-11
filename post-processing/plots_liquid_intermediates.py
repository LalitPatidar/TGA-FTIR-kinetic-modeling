#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os

time = np.asarray([])
temperature = np.asarray([])
mass_fraction = np.asarray([])
t_HONO = np.asarray([])
c_HONO = np.asarray([])
HCOOH = np.asarray([])
N2 = np.asarray([])
ONNO2 = np.asarray([])
ONTNTA = np.asarray([])
CH2NNO2 = np.asarray([])
INT86a = np.asarray([])
INT249a = np.asarray([])
CH2NH = np.asarray([])

with open('./mass_fractions_liquid.txt','r') as File:
    lines = File.readlines()
    
    header = lines[0].split()
    
    it_HONO = header.index('t_HONO')
    ic_HONO = header.index('c_HONO')
    iHCOOH = header.index('HCOOH')
    iN2 = header.index('N2')
    iONNO2 = header.index('ONNO2')
    iONTNTA = header.index('ONTNTA')
    iCH2NNO2 = header.index('CH2NNO2')
    iINT86a = header.index('INT86a')
    iINT249a = header.index('INT249a')
    iCH2NH = header.index('CH2NH')
    
    new_lines = lines[1::2]

    for line in new_lines:
        time = np.append(time,float(line.split()[0]))
        temperature = np.append(temperature,float(line.split()[1])-273.0)
        mass_fraction = np.append(mass_fraction,100.0*float(line.split()[2]))
        mj = float(line.split()[2])*1.492
        t_HONO = np.append(t_HONO,mj*float(line.split()[it_HONO]))
        c_HONO = np.append(c_HONO,mj*float(line.split()[ic_HONO]))
        HCOOH = np.append(HCOOH,mj*float(line.split()[iHCOOH]))
        N2 = np.append(N2,mj*float(line.split()[iN2]))
        ONNO2 = np.append(ONNO2,mj*float(line.split()[iONNO2]))
        ONTNTA = np.append(ONTNTA,mj*float(line.split()[iONTNTA]))
        CH2NNO2 = np.append(CH2NNO2,mj*float(line.split()[iCH2NNO2]))
        INT86a = np.append(INT86a,mj*float(line.split()[iINT86a]))
        INT249a = np.append(INT249a,mj*float(line.split()[iINT249a]))
        CH2NH = np.append(CH2NH,mj*float(line.split()[iCH2NH]))

T_start = 275.0
T_end = 300.0

fig, ax1 = plt.subplots()
"""
ax1.plot(temperature,HMX,label='HMX',marker='*',color='black')
ax1.set_ylabel('Mass of HMX in liquid (mg)')
ax1.legend(loc='upper left',bbox_to_anchor=(0.0,1.15),ncol=1)
#ax1.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax1.set_xlim(T_start,T_end)
ax1.set_xlabel('Temperature($^o$C)')

ax2 = ax1.twinx()
"""
#ax1.plot(temperature,t_HONO,label='tHONO',marker=">")
#ax1.plot(temperature,c_HONO,label='cHONO',marker="p")
#ax1.plot(temperature,HCOOH,label='HCOOH',marker="<")
#ax1.plot(temperature,N2,label='N2',marker="D")
ax1.plot(temperature,ONNO2,label='ONNO2',marker="^")
ax1.plot(temperature,ONTNTA,label='ONTNTA',marker="d")
ax1.plot(temperature,CH2NNO2,label='CH2NNO2',marker="o")
ax1.plot(temperature,INT86a,label='INT86a',marker="s")
ax1.plot(temperature,INT249a,label='INT249a',marker="d")
ax1.plot(temperature,CH2NH,label='CH2NH',marker="d")

#ax2.ticklabel_format(style='sci',axis='y',scilimits=(-1,1))
ax1.legend(loc='upper center',bbox_to_anchor=(0.5,1.15),ncol=3)
ax1.set_ylabel('Mass of intermediate species in liquid (mg)')

plt.savefig('Liquids_intermediate.pdf')
