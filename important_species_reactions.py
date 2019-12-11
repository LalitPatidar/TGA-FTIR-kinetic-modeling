import os
import numpy as np

    
path_input = './input'

with open(os.path.join(path_input,'nitramine-liquid-phase-mechanism.txt'),'r') as File:
    lines = File.readlines()
    
    elementsList = []
    copy = False
    
    for line in lines:
        if line.strip() == 'ELEMENTS':
            copy = True
        if line.strip() == 'END':
            copy = False
        elif copy and line.strip() != 'ELEMENTS':
            elementsList.append(line.strip())
    elementsList.remove('Ar')
    

    thermoList = np.asarray([])
    copy = False
    
    for line in lines:
        if line.strip() == 'THERMO':
            copy = True
        if line.strip() == 'END':
            copy = False
        elif copy and line.strip() != 'THERMO':
            thermoList = np.append(thermoList,line)
        
    speciesList = [thermoList[i] for i in range(0, len(thermoList), 4)]

    reactionsList = []
    copy = False

    for line in lines:
        if line.strip() == 'REACTIONS':
            copy = True
        if line.strip() == 'END':
            copy = False
        elif copy and line.strip() != 'REACTIONS':
            reactionsList.append(line.strip().split()[0])
           
speciesData = []
for item in speciesList:
    line1 = item.split()
    name = line1[0]
    
    compositionDict = {}
    compositionDict[line1[1]] = int(line1[2][0])
    compositionDict[line1[2][1]] = int(line1[3][0])
    compositionDict[line1[3][1]] = int(line1[4][0])
    compositionDict[line1[4][1]] = int(line1[5][0])
    
    speciesData.append([name,compositionDict])
    
reactionsData = []
for item in reactionsList:
    reactants_list = item.split('=')[0].split('+')
    products_list = item.split('=')[1].split('+')
    
    reactionsData.append([item,reactants_list,products_list])
#--------------------------------------------------------------------------------

path_input = './output_TGA'

with open(os.path.join(path_input,'rates_of_production.txt'),'r') as File:
    lines = File.readlines()
      
    header = lines[0].split()
    
    array1 = np.asarray([float(y) for y in lines[1].split()])
    
    for line in lines[2:]:
        #temp = np.asarray([float(y) for y in line.split()])
        temp = np.asarray([abs(float(y)) for y in line.split()])
        array1 = np.vstack((array1,temp))
       

dt = 2.0
Vt0 = 1.492e-3/1.8 # volume in cm3 = mass(g)/density(g/cm3)

for i in range(3,array1.shape[1]):
    array1[:,i] = dt*Vt0*np.multiply(array1[:,2],array1[:,i])
        
W_np = np.sum(array1,axis=0)
W_list = [[W_np[i], header[i]] for i in range(3,len(W_np))]

W_final = sorted(W_list,key=lambda x:x[0], reverse=True)

with open(os.path.join(path_input,'important_species.txt'),'w') as File:
    for item in W_final:
        File.writelines('%50s %15.2E\n' %(item[1],item[0]))
#----------------------------------------------------------------------------------------    
    
with open(os.path.join(path_input,'rates_of_progress.txt'),'r') as File:
    lines = File.readlines()
      
    header = lines[0].split()
    
    array1 = np.asarray([float(y) for y in lines[1].split()])
    
    for line in lines[2:]:
        #temp = np.asarray([float(y) for y in line.split()])
        temp = np.asarray([abs(float(y)) for y in line.split()])
        array1 = np.vstack((array1,temp))
       

dt = 2.0
Vt0 = 1.492e-3/1.8 # volume in cm3 = mass(g)/density(g/cm3)

for i in range(3,array1.shape[1]):
    array1[:,i] = dt*Vt0*np.multiply(array1[:,2],array1[:,i])
        
q_np = np.sum(array1,axis=0)
q_list = [[q_np[i], reactionsData[i-3][0]] for i in range(3,len(q_np))]

q_final = sorted(q_list,key=lambda x:x[0], reverse=True)

with open(os.path.join(path_input,'important_reactions.txt'),'w') as File:
    for item in q_final:
        File.writelines('%50s %15.2E\n' %(item[1],item[0]))
        
