import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator


temperature = np.asarray([])
H2CO = np.asarray([])
CO = np.asarray([])
HCN = np.asarray([])
N2O = np.asarray([])
NO = np.asarray([])
NO2 = np.asarray([])
H2O = np.asarray([])
CO2 = np.asarray([])

Path = './partial_pressures'
file_out = 'HMX_15KPM_1492ug'

start = 1
step = 20

heating_rate = 15.0
OPUS_start_temp = 250.0
time_resolution = 0.144
#T_offset = heating_rate*16.95/60
T_offset = 0.0

with open(Path + '/H2CO.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step] 
    for line in lines:
        temperature = np.append(temperature,time_resolution*float(line.split()[0])*(heating_rate/60.0) + OPUS_start_temp - T_offset)
        H2CO = np.append(H2CO,float(line.split()[2]))
        
with open(Path + '/CO.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        CO = np.append(CO,float(line.split()[2]))
        
with open(Path + '/HCN.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
    for line in lines:
        HCN = np.append(HCN,float(line.split()[2]))
        
with open(Path + '/N2O.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
  
    for line in lines:
        N2O = np.append(N2O,float(line.split()[2]))
        
with open(Path + '/NO.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
  
    for line in lines:
        NO = np.append(NO,float(line.split()[2]))
        
with open(Path + '/NO2.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
  
    for line in lines:
        NO2 = np.append(NO2,float(line.split()[2]))
        
with open(Path + '/H2O.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
  
    for line in lines:
        H2O = np.append(H2O,float(line.split()[2]))
        
with open(Path + '/CO2.txt','r') as File:
    lines = File.readlines()
    lines = lines[start::step]
  
    for line in lines:
        CO2 = np.append(CO2,float(line.split()[2]))

fig, ax1 = plt.subplots()
ax1.xaxis.set_minor_locator(AutoMinorLocator())

ax1.plot(temperature,H2CO,label='CH$_2$O',marker=">")
ax1.plot(temperature,CO,label='CO',marker="p")
ax1.plot(temperature,HCN,label='HCN',marker="<")
ax1.plot(temperature,N2O,label='N$_2$O',marker="D")
ax1.plot(temperature,NO,label='NO',marker="^")
ax1.plot(temperature,NO2,label='NO$_2$',marker="o")
ax1.plot(temperature,H2O,label='H$_2$O',marker="s")
ax1.plot(temperature,CO2,label='CO$_2$',marker="d")

ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
ax1.set_xlim(275,300)
ax1.set_xlabel('Temperature($^o$C)')
ax1.set_ylabel('Mole fraction (%)')
ax1.legend(loc='upper center',bbox_to_anchor=(0.5,1.15),ncol=4)

ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(direction='in',which='both')

textstr = '\n'.join(['','Experimental'])
ax1.text(0.75,0.92,textstr,transform=ax1.transAxes)

plt.savefig('Gases-fig9'+file_out+'.pdf')
