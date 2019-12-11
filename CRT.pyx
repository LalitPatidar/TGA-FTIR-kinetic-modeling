# cython: language_level=3, boundscheck=False

import os
import numpy as np
cimport numpy as np
np.import_array()
import scipy.integrate
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Phase.h":
    cdef cppclass Phase:
        Phase(string,string) except +
        double getT()
        void setT(double)
        double getP()
        void setP(double)
        double get_density()
        void set_density(double)
        double get_viscosity()
        void set_viscosity(double)
                
        int get_n_species()
        int get_n_reactions()
        string species_name(int)
        string reaction_name(int)
        int species_index(string)
        vector[string] get_reactants_list(int)
        vector[string] get_products_list(int)
        
        void setX(vector[double])
        vector[double] getX()
        void setY(vector[double])
        vector[double] getY()
        
        vector[double] get_molecular_weights()
        double get_mean_molecular_weight()
        vector[double] get_concentrations()
        
        vector[double] get_forward_rate_constants()
        vector[double] get_backward_rate_constants()
        
        vector[double] get_kdiff_f()
        vector[double] get_kdiff_b()
        
        vector[vector[double]] get_kdiff()
        
        vector[double] get_net_forward_rate_constants()
        vector[double] get_net_backward_rate_constants()

        vector[double] get_Eaf()
        vector[double] get_Eab()
        
        vector[double] get_net_rates_of_production()
        vector[double] get_net_rates_of_progress()              
        
cdef class PyPhase:
    cdef Phase* liquidphase
        
    cdef np.ndarray _kdiff_f, _kdiff_b

    def __cinit__(self,file1,file2):
        self.liquidphase = new Phase(file1,file2)
        
        
    @property    
    def T(self):
        return self.liquidphase.getT()
        
    @T.setter    
    def T(self,double value):
        self.liquidphase.setT(value)
        self.liquidphase.set_viscosity(self.viscosity)
        self.liquidphase.set_density(self.density)
        
    @property    
    def P(self):
        return self.liquidphase.getP()
        
    @P.setter    
    def P(self,double value):
        self.liquidphase.setP(value)
        
    @property    
    def density(self):
        T = [298.15, 500.0, 550.0, 600.0, 650.0,700.0, 750.0, 800.0]
        rho = [1.8, 1.8, 1.6509, 1.6144, 1.5869, 1.5545, 1.5201, 1.4882]
        # density in g/cm3
        return np.interp(self.T,T,rho)
        #return 1.8
        
    @density.setter    
    def density(self,double value):
        self.liquidphase.set_density(value)
        
    @property    
    def viscosity(self):
        #return self.liquidphase.get_viscosity()
        T = [298.15,550.0,600.0,650.0,700.0,750.0,800.0]
        eta = [0.45, 0.45,0.12,0.04,0.0220,0.01,0.0055]
        # viscosity in g/cm-s
        return 10.0*np.interp(self.T,T,eta)
        #return 0.12
        
    @viscosity.setter
    def viscosity(self,double value):
        self.liquidphase.set_viscosity(value)
                        
    @property
    def n_species(self):
        return self.liquidphase.get_n_species()
        
    @property
    def n_reactions(self):
        return self.liquidphase.get_n_reactions()
        
    def species_name(self, int k):
        return self.liquidphase.species_name(k).decode('ASCII')
        
    def reaction_name(self, int k):
        return self.liquidphase.reaction_name(k).decode('ASCII')
        
    def species_index(self, name):
        return self.liquidphase.species_index(name.encode('utf-8'))
    
    @property
    def Eaf(self):
        return self.liquidphase.get_Eaf()
        
    @property
    def Eab(self):
        return self.liquidphase.get_Eab()
    
    @property
    def X(self):
        return self.liquidphase.getX()
        
    @X.setter
    def X(self, np.ndarray values):
        self.liquidphase.setX(values)
        
    @property
    def Y(self):
        return self.liquidphase.getY()
        
    @Y.setter
    def Y(self, np.ndarray values):
        self.liquidphase.setY(values)
        
    @property
    def molecular_weights(self):
        return self.liquidphase.get_molecular_weights()
        
    @property
    def mean_molecular_weight(self):
        return self.liquidphase.get_mean_molecular_weight()
        
    @property
    def concentrations(self):
        return self.liquidphase.get_concentrations()    
            
    @property
    def forward_rate_constants(self):
        return self.liquidphase.get_forward_rate_constants()
        
    @property
    def backward_rate_constants(self):
        return self.liquidphase.get_backward_rate_constants()

    @property
    def kdiff_f(self):
        return self.liquidphase.get_kdiff_f()
        
    @property
    def kdiff_b(self):
        return self.liquidphase.get_kdiff_b()
        
    @property
    def kdiff(self):
        return np.asarray(self.liquidphase.get_kdiff())
        
    @property
    def net_forward_rate_constants(self):
        return self.liquidphase.get_net_forward_rate_constants()
        
    @property
    def net_backward_rate_constants(self):
        return self.liquidphase.get_net_backward_rate_constants()
        
    def reactants_list(self,int k):
        return [item.decode('utf-8') for item in list(self.liquidphase.get_reactants_list(k))]
        
    def products_list(self,int k):
        return [item.decode('utf-8') for item in list(self.liquidphase.get_products_list(k))]
        
    @property
    def net_rates_of_production(self):
        return self.liquidphase.get_net_rates_of_production()
        
    @property
    def net_rates_of_progress(self):
        return self.liquidphase.get_net_rates_of_progress()
    
####################################################################################################################################
cdef class ReactorOde(object):
    cdef object liquid
    cdef float Tinit, Tend
    cdef dict evaporating_species, gas_cell_parameters
    cdef float Pg, Vg, Tg, mdot_purge
    cdef list global_AnE
    cdef int evaporation_global
    
    def __init__(self, liquid, Tinit, Tend, evaporating_species, gas_cell_parameters, global_AnE, evaporation_global):
        self.liquid = liquid
        self.Tinit = Tinit
        self.Tend = Tend
        self.evaporating_species = evaporating_species
        self.Pg = gas_cell_parameters["Gas-cell-Pressure(Pa)"]
        self.Vg = gas_cell_parameters["Gas-cell-Volume(m3)"]
        self.Tg = gas_cell_parameters["Gas-cell-Temperature(K)"]
        self.mdot_purge = gas_cell_parameters["Purge-gas-flow-rate(g/s)"]
        self.global_AnE = global_AnE
        self.evaporation_global = evaporation_global

        
    def __call__(self, float t, np.ndarray[np.float64_t,ndim=1] y):
        """the ODE function, y' = f(t,y) """
        
        cdef int n_species = self.liquid.n_species

        cdef float Pg = self.Pg
        cdef float Vg = self.Vg
        cdef float Tg =  self.Tg
        cdef float mdot_purge = self.mdot_purge
        cdef float MW_purge = 28.014 #g/mol
        cdef float R = 8.314
        
        # Default evaporation parameters for all species
        cdef np.ndarray A = np.asarray([0.0 for i in range(n_species)])
        cdef np.ndarray n = np.asarray([0.0 for i in range(n_species)])
        cdef np.ndarray Ea = np.asarray([0.0 for i in range(n_species)])
        #cdef np.ndarray[np.float64_t,ndim=1] k_lg = np.asarray([0.0 for i in range(n_species)])
        
        #"""
        # Assisgn evaporation-parameters
        for name, local_AnEa in self.evaporating_species.items():
            if self.evaporation_global:
                A_factor = self.global_AnE[0]
                n_factor = self.global_AnE[1]
                Ea_factor = self.global_AnE[2]
            else:
                A_factor = local_AnEa[0]
                n_factor = local_AnEa[1]
                Ea_factor = local_AnEa[2]
                
            A[self.liquid.species_index(name)] = A_factor
            n[self.liquid.species_index(name)] = n_factor
            Ea[self.liquid.species_index(name)] = Ea_factor

        #"""
            
        self.liquid.T = y[0]

        y = np.asarray([max([y[i],0.0]) for i in range(len(y))])

        y[2:2+n_species] = y[2:2+n_species]/sum(y[2:2+n_species])
        self.liquid.Y = y[2:2+n_species]
        
        y[2+n_species:2+n_species+n_species] = y[2+n_species:2+n_species+n_species]/sum(y[2+n_species:2+n_species+n_species])
        
        # Note different unit of R here
        k_lg = np.multiply(np.power(self.liquid.T,n),np.multiply(A,np.exp(-Ea/1.9872/self.liquid.T)))
        
        cdef np.ndarray[np.float64_t,ndim=1] KY = np.multiply(k_lg,y[2:2+n_species])
        
        cdef np.ndarray[np.float64_t,ndim=1] wdot = np.asarray(self.liquid.net_rates_of_production)
        
        #cdef float dTdt = self.Heating_rate/60.0
        cdef float dTdt = (self.Tend - self.Tinit)*22.2*np.exp(-22.2*t)
            
        #cdef float dYcdt = -y[1]*np.sum(KY)
        cdef float dMcdt = -y[1]*np.sum(KY)
        
        cdef np.ndarray[np.float64_t,ndim=1] dYjcdt = np.multiply(wdot,self.liquid.molecular_weights)/self.liquid.density - KY \
                + np.multiply(y[2:2+n_species],np.sum(KY))
                
        cdef np.ndarray[np.float64_t,ndim=1] dXgdt = (R*Tg/Pg/Vg)*\
                (y[1]*np.divide(KY[:n_species-1],self.liquid.molecular_weights[:n_species-1]) \
                 -np.multiply((np.sum(np.divide(y[1]*KY,self.liquid.molecular_weights)) + mdot_purge/MW_purge),y[2+n_species:2+n_species+n_species-1]))
                
        cdef float dXpdt = (R*Tg/Pg/Vg)*\
                (mdot_purge/MW_purge + y[1]*np.divide(KY[n_species-1],self.liquid.molecular_weights[n_species-1])\
                 -np.multiply((np.sum(np.divide(y[1]*KY,self.liquid.molecular_weights)) + mdot_purge/MW_purge),y[2+n_species+n_species-1]))
        
        return np.hstack((dTdt,dMcdt, dYjcdt, dXgdt, dXpdt))

####################################################################################################################################

cpdef main(parameters):
    #------------------------
    # Unpack input parameters
    #------------------------
    path_input = parameters["path_input"]
    path_output = parameters["path_output"]
    
    Tinit = parameters["Tinit"] + 273.15
    Tend = parameters["Tend"] + 273.15
    
    t_end = parameters["t_end"]

    dt = parameters["dt"]
    Mc0 = 0.001*parameters["Mc0"]
    
    file_liquid_phase_mechanism = './input/' + parameters["liquid-phase-mechanism"]
    file_species_radius = './input/' + parameters["species_radius"]
    file_evaporation_parameters = parameters["evaporation_parameters"]
    evaporation_global = parameters["evaporation_global"] 
    
    Reactant = parameters["Reactant"]
    
    gas_cell_parameters = {k:parameters[k] for k in ("Gas-cell-Pressure(Pa)","Gas-cell-Volume(m3)","Gas-cell-Temperature(K)","Purge-gas-flow-rate(g/s)") if k in parameters}   

    liquid = PyPhase(file_liquid_phase_mechanism.encode('utf-8'),file_species_radius.encode('utf-8'))
    
    #---------------------------------------------------------------------
    # Initialize with temperature and density 
    # temperature dependent viscosity is implement via temperature setter
    # similarly temperature dependent density could be implemented
    #---------------------------------------------------------------------
    liquid.T = Tinit
    cdef int n_species = liquid.n_species
    cdef int n_reactions = liquid.n_reactions
    cdef int i = 0
    liquid.X = np.asarray([1.0 if liquid.species_name(i)==Reactant else 0.0 for i in range(n_species)])

    #Yc0 = 1.0
    Xg = np.zeros(n_species-1)
    Xp = 1.0

    #y0 = np.hstack((liquid.T, Yc0, liquid.Y, Xg, Xp))
    y0 = np.hstack((liquid.T, Mc0, liquid.Y, Xg, Xp))
    
    # Read evaporation-parameters
    cdef dict evaporating_species = {}
    with open(os.path.join(path_input,file_evaporation_parameters),'r') as File:
        lines = File.readlines()
            
        for line in lines:
            if line.startswith('!'):
                pass
            else:
                name = line.split()[0]
                evaporating_species[name] = [float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
    
    # Default evaporation parameters for all species
    cdef np.ndarray A = np.asarray([0.0 for i in range(n_species)])
    cdef np.ndarray n = np.asarray([0.0 for i in range(n_species)])
    cdef np.ndarray Ea = np.asarray([0.0 for i in range(n_species)])
    cdef np.ndarray klg = np.asarray([0.0 for i in range(n_species)])
        
    # Assisgn evaporation-parameters for evaporating species
    global_AnE = parameters["global_AnE"]
    for name, local_AnEa in evaporating_species.items():
        if evaporation_global:
            A_factor = global_AnE[0]
            n_factor = global_AnE[1]
            Ea_factor = global_AnE[2]
        else:
            A_factor = local_AnEa[0]
            n_factor = local_AnEa[1]
            Ea_factor = local_AnEa[2]
            
        A[liquid.species_index(name)] = A_factor
        n[liquid.species_index(name)] = n_factor
        Ea[liquid.species_index(name)] = Ea_factor   

    # Set up objects representing the ODE and the solver
    ode = ReactorOde(liquid,Tinit,Tend, evaporating_species, gas_cell_parameters, global_AnE, evaporation_global)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-12,rtol=1e-6,nsteps=10000)
    solver.set_initial_value(y0, 0.0)
    
    # Integrate the equations, keeping T(t) and Y(k,t)
    
    mass_fractions_liquid = []
    mole_fractions_gas = []

    iReactant = liquid.species_index(Reactant)
    #------------------------------------------------------
    header = '%15s %15s %15s' % ('Time(s)','Temp(K)', 'Mc(t)/Mc(o)')
    
    for i in range(n_species):
        header = header + '%15s' %(liquid.species_name(i))
    header = header + '\n'
    
    with open(os.path.join(path_output,'mass_fractions_liquid.txt'),'w') as File:
        File.writelines(header)
    #------------------------------------------------------    
    header = '%15s %15s %15s' % ('Time(s)','Temp(K)', 'Mc(t)/Mc(o)') 
    
    for i in range(n_species):
        header = header + '%15s' %(liquid.species_name(i))
    header = header + '\n'
          
    with open(os.path.join(path_output,'mole_fractions_gas.txt'),'w') as File:
        File.writelines(header)
    #------------------------------------------------------    
    """
    header = '%15s' %('Time(s)')
    for i in range(n_reactions):
        header = header + '%15d' %(i+1)
    header = header + '\n'
     
    with open(path2output+'rates_of_progress.txt','w') as File:
        #File.writelines(header)
        File.writelines("Top 10 sensitive reactions\n")
    global_q = np.zeros(n_reactions)
    #------------------------------------------------------
    with open(path2output+'rates_of_production.txt','w') as File:
        #File.writelines(header)
        File.writelines("Top 10 sensitive species\n")
    #------------------------------------------------------
    with open(path2output+'rates_of_evaporating_species.txt','w') as File:
        #File.writelines(header)
        File.writelines("Rates of chemical production and evaporation of species\n")
    #------------------------------------------------------
    """
    while solver.successful() and solver.t < t_end:    
        liquid.T = solver.y[0]
        #Yc = solver.y[1]
        Mc = solver.y[1]
        Yc = Mc/Mc0
        
        Yjc = np.asarray([max([solver.y[i],0.0]) for i in range(len(solver.y))])
        
        liquid.Y = Yjc[2:2+n_species]/sum(Yjc[2:2+n_species])
        
        Xjg = Yjc[2+n_species:2+n_species+n_species]/sum(Yjc[2+n_species:2+n_species+n_species])
        
        if solver.t%0.01<dt:
            line_liquid = '%15.3f %15.3f %15.3E' %(solver.t, solver.y[0], Yc)
            line_gas = '%15.3f %15.3f %15.3E' % (solver.t, solver.y[0], Yc)
            line_reactions = '%15.3f' % (solver.y[0]-273.15)
            
            for i in range(n_species):
                line_liquid = line_liquid + '%15.3E' %(liquid.Y[i])
                line_gas = line_gas + '%15.3E' %(Xjg[i])
                
            line_liquid = line_liquid + '\n'
            line_gas = line_gas + '\n'
            
            for i in range(n_reactions):
                #global_q[i] = global_q[i] + abs(liquid.net_rates_of_progress[i])
                #line_reactions = line_reactions + '%50.3E' %(global_q[i])
                line_reactions = line_reactions + '%15.3E' %(liquid.net_rates_of_progress[i])
                
            line_reactions = line_reactions + '\n'

            print ('%15s %15.5f %15s %15.5f %15s %15.3E %15s %15.3E' %('time(s) = ', solver.t,'temperature(oC) = ', liquid.T - 273.15, 'Y_HMX = ', liquid.Y[iReactant], 'Yc = ', Yc))
            
            with open(os.path.join(path_output,'mass_fractions_liquid.txt'),'a') as File:
                File.writelines(line_liquid)
            
            """    
            with open(os.path.join(path_output,'rates_of_progress.txt'),'a') as File:
                #File.writelines(line_reactions)
                File.writelines('%15.3f\n' % solver.t)
                for indx in np.absolute(liquid.net_rates_of_progress).argsort()[-30:][::-1]:
                    File.writelines('%50s %15.3E' %(liquid.reaction_name(indx),liquid.net_rates_of_progress[indx]))
                    File.writelines("\n")
                File.writelines("-----------------------------------------------------------------\n")
            
                
            with open(os.path.join(path_output,'rates_of_production.txt'),'a') as File:
                File.writelines('%15.3f\n' % solver.t)
                for indx in np.absolute(liquid.net_rates_of_production).argsort()[-10:][::-1]:
                    File.writelines('%20s %15.3E' %(liquid.species_name(indx),liquid.net_rates_of_production[indx]))
                    File.writelines("\n")
                File.writelines("-----------------------------------------------------------------\n")           
            
            with open(os.path.join(path_output,'rates_of_evaporating_species.txt'),'a') as File:
                File.writelines('%-10s %-5.3f %15s %15s %15s %15s\n' % ('Temp(C)=',solver.y[0]-273.0, 'Y_liquid', 'wdot_chem', 'wdot_evap', 'klg'))
                for name in evaporating_species:
                    indx = liquid.species_index(name)
                    klg = np.multiply(np.power(liquid.T,n),np.multiply(A,np.exp(-Ea/1.9872/liquid.T)))
                    ndot_evap = (liquid.density*liquid.Y[indx]/liquid.molecular_weights[indx])*klg[indx]
                    File.writelines('%20s %15.3E %15.3E %15.3E %15.3E' %(liquid.species_name(indx), liquid.Y[indx], liquid.net_rates_of_production[indx], ndot_evap, klg[indx]))
                    File.writelines("\n")
                File.writelines("------------------------------------------------------------------------------------------\n")
            """       
            with open(os.path.join(path_output,'mole_fractions_gas.txt'),'a') as File:
                File.writelines(line_gas)
            
        solver.integrate(solver.t + dt)
    """    
    with open(os.path.join(path_output,'rates_of_progress.txt'),'a') as File:
        #File.writelines(line_reactions)
        File.writelines('%15.3f\n' % solver.t)
        for indx in np.absolute(liquid.net_rates_of_progress).argsort()[-30:][::-1]:
            File.writelines('%50s %15.3E' %(liquid.reaction_name(indx),liquid.net_rates_of_progress[indx]))
            File.writelines("\n")
        File.writelines("-----------------------------------------------------------------\n")
                
    with open(os.path.join(path_output,'rates_of_production.txt'),'a') as File:
        File.writelines('%15.3f\n' % solver.t)
        for indx in np.absolute(liquid.net_rates_of_production).argsort()[-10:][::-1]:
            File.writelines('%20s %15.3E' %(liquid.species_name(indx),liquid.net_rates_of_production[indx]))
            File.writelines("\n")
        File.writelines("-----------------------------------------------------------------\n")
    
    with open(os.path.join(path_output,'rates_of_evaporating_species.txt'),'a') as File:
        File.writelines('%-10s %-5.3f %15s %15s %15s %15s\n' % ('Temp(C)=',solver.y[0]-273.0, 'Y_liquid', 'wdot_chem', 'wdot_evap', 'klg'))
        for name in evaporating_species:
            indx = liquid.species_index(name)
            klg = np.multiply(np.power(liquid.T,n),np.multiply(A,np.exp(-Ea/1.9872/liquid.T)))
            ndot_evap = (liquid.density*liquid.Y[indx]/liquid.molecular_weights[indx])*klg[indx]
            File.writelines('%20s %15.3E %15.3E %15.3E %15.3E' %(liquid.species_name(indx), liquid.Y[indx], liquid.net_rates_of_production[indx], ndot_evap, klg[indx]))
            File.writelines("\n")
        File.writelines("------------------------------------------------------------------------------------------\n")
    """
    print('No. of species = %d'%(liquid.n_species))
    print('No. of reactions = %d'%(liquid.n_reactions))
    
