import matplotlib.pyplot as plt
import numpy as np

time = np.asarray([])
temperature = np.asarray([])
CH2O = np.asarray([])
CO = np.asarray([])
HCN = np.asarray([])
N2O = np.asarray([])
NO = np.asarray([])
NO2 = np.asarray([])
H2O = np.asarray([])
CO2 = np.asarray([])

with open('./output_CRT/mole_fractions_gas.txt','r') as File:
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
    
    new_lines = lines[1:70:50]
    temp = lines[70::100]
    
    new_lines.extend(temp)
    
    for line in new_lines:
        time = np.append(time,float(line.split()[0]))
        temperature = np.append(temperature,float(line.split()[1])-273.0)
        CH2O = np.append(CH2O,float(line.split()[iCH2O]))
        CO = np.append(CO,float(line.split()[iCO]))
        HCN = np.append(HCN,float(line.split()[iHCN]))
        N2O = np.append(N2O,float(line.split()[iN2O]))
        NO = np.append(NO,float(line.split()[iNO]))
        NO2 = np.append(NO2,float(line.split()[iNO2]))
        H2O = np.append(H2O,float(line.split()[iH2O]))
        CO2 = np.append(CO2,float(line.split()[iCO2]))
        

t_end = 30.0

plt.figure(1)
plt.plot(time,CH2O/2.5e-3,label='CH2O',marker=">")
plt.plot(time,CO/2.5e-3,label='CO',marker="p")
plt.plot(time,HCN/2.5e-3,label='HCN',marker="<")
plt.plot(time,N2O/2.5e-3,label='N2O',marker="D")
plt.plot(time,NO/2.5e-3,label='NO',marker="^")
plt.plot(time,NO2/2.5e-3,label='NO2',marker="o")
plt.plot(time,H2O/2.5e-3,label='H2O',marker="s")
plt.plot(time,CO2/2.5e-3,label='CO2',marker="d")
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.xlim(0.0,t_end)
plt.xlabel('Time(s)')
plt.ylabel('Relative Mole fraction')
plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.15),ncol=4)
plt.savefig('./output_CRT/Gases.pdf')

plt.figure(2)
plt.plot(time,temperature)
plt.ylabel('Temperature($^o$C)')
plt.xlabel('time (s)')
plt.xlim(0.0,1.0)
plt.savefig('./output_CRT/Temperature.pdf')
