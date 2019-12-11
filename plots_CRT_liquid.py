import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


with open('./output_CRT/mass_fractions_liquid.txt','r') as File:
    lines = File.readlines()
    
    header = lines[0].split()
    
    array1 = np.asarray([float(y) for y in lines[1].split()])
    
    for line in lines[2:]:
        temp = np.asarray([float(y) for y in line.split()])
        array1 = np.vstack((array1,temp))

with PdfPages('./output_CRT/check_plots.pdf') as pdf:
    for i in range(3,array1.shape[1]):
        fig, ax1 = plt.subplots()
        ax1.plot(array1[:,0],1.492*np.multiply(array1[:,2],array1[:,i]),label=header[i])

        ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Mass in liquid (mg)')
        ax1.legend(loc='upper center',bbox_to_anchor=(0.5,1.15))
        #plt.savefig('./plots_check/'+ header[i] +'.pdf')
        pdf.savefig()
        plt.close()
    
