import os
import time
import numpy as np

import CRT

start = time.time()

path_input = './input'
path_output = './output_CRT'

#---------------- Automatic mechanism reader-----------------------------------------------------------
#--------------- (to be removed once kinetics and evaporation are optimized)---------------------------
#------------------------------------------------------------------------------------------------------
with open(os.path.join(path_input,'nitramine-liquid-phase-mechanism-input.txt'),'r') as File:
    lines = File.readlines()
    
    flag1 = 0
    species_list = []
    reactions_list = []
    
    for line in lines:
        if line.startswith('REACTIONS'):
            flag1 = 1

        if (flag1 == 1 and (line.strip() not in ['!','REACTIONS','END'])) and not line.startswith('!'):
            reactions_list.append(line)
            reactants,products = line.split(' ')[0].split('=')
            species = reactants.split('+') + products.split('+')
            species_list.extend(species)

seen = set()
seen_add = seen.add
species_list = [x for x in species_list if not(x in seen or seen_add(x))]

with open(os.path.join(path_input,'nitramine-liquid-phase-mechanism.txt'),'w') as File2:
    File2.writelines('ELEMENTS\nC\nH\nN\nO\nAr\nEND\n')
    File2.writelines('SPECIES\n')
    for species in species_list:
        File2.writelines(species + '\n')
    File2.writelines('N2\n')
    File2.writelines('END\n')
    File2.writelines('THERMO\n')

with open(os.path.join(path_input,'thermo-data.txt'),'r') as File:
    lines = File.readlines()
    
    
    for species in species_list:    
        for i in range(len(lines)):
            if lines[i].split()[0] == species and lines[i].split()[0] != 'N2':
                with open(os.path.join(path_input,'nitramine-liquid-phase-mechanism.txt'),'a') as File2:
                    File2.writelines(lines[i:i+4])
                    
            if lines[i].split()[0] == 'N2':
                N2_lines = lines[i:i+4]
                    
    with open(os.path.join(path_input,'nitramine-liquid-phase-mechanism.txt'),'a') as File2:
        File2.writelines(N2_lines)
        
with open(os.path.join(path_input,'nitramine-liquid-phase-mechanism.txt'),'a') as File2:
    File2.writelines('END\n')
    File2.writelines('REACTIONS\n')
    File2.writelines(reactions_list)
    File2.writelines('END')
    
#------------------------------------------------------------------------------

#---------------- Input parameters---------------------------------------------
parameters = {
                "Reactant":'HMX',
                "liquid-phase-mechanism":'nitramine-liquid-phase-mechanism.txt',
                "species_radius":'log-file-data-minima.txt',
                "evaporation_parameters":'evaporation-parameters.txt',
                
                "evaporation_global":1,
                "global_AnE":[6.38e18,0,51.3e3],
                
                "Tinit":250.00,
                "Tend":300.00,
                "t_end":50.0,
                "dt":0.0001,
                "Mc0":1.0,
                
                "Gas-cell-Pressure(Pa)":101325.0,
                "Gas-cell-Volume(m3)":8.7e-6,
                "Gas-cell-Temperature(K)":473.15,
                "Purge-gas-flow-rate(g/s)":70.0*(1.13*0.001/60.0),
                
                "path_input":path_input,
                "path_output":path_output
}
#------------------------------------------------------------------------------
    
CRT.main(parameters)

end = time.time()

print('Wall time = %f'%(end-start))
