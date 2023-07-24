import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sp 
from Reactor_MODEL import Reactor
from UTSG_MODEL import UTSG
from reactor_pump import pump
import threading 

class Pressurizer():

    def __init__(self,area:float,water_enthalpy:float,latent_heat_of_vaporization:float,
                 water_enthalpy_spray_nozzle:float,heat_conversion_factor:float,
                 effective_length_of_pressurizer:float,water_level:float,
                 mass_steam:float,mass_water:float,pressure_at_pressurizer:float,
                 condensation_flow_rate:float,spray_flow_rate:float,in_and_out_flow_rate:float,
                 water_vol_in_pressurizer:float,specific_vol_water:float,specific_vol_steam:float,
                 steam_density:float,water_density:float,pressure_cont:list):
        
        #design parameters
        self.area=area
        self.L=effective_length_of_pressurizer
        self.J=heat_conversion_factor

        #constants 
        self.const=np.array(pressure_cont)

        """ const[0]=d meu_f/dpr
            const[1]=d meu_g/dpr
            const[2]=d meu_fg/dpr
            const[3]=d rou_s/dpr
            const[4]=d rou_w/dpr
            const[5]=d h_f/dpr
            const[6]=d h_fg/dpr 
        """

        self.h_f=water_enthalpy
        self.h_fg=latent_heat_of_vaporization
        self.h_sp=water_enthalpy_spray_nozzle
        self.rou_w=water_density
        self.rou_s=steam_density
        self.meu_f=specific_vol_water
        self.meu_g=specific_vol_steam
        self.W_sp=spray_flow_rate
        self.W_sr=in_and_out_flow_rate
        

        #variables
        self.M_s=mass_steam
        self.M_w=mass_water
        self.L_w=water_level
        self.V_w=water_vol_in_pressurizer
        self.P_pr=pressure_at_pressurizer
        self.W_co=condensation_flow_rate

        #calculating some constant parameters 

        self.C1=self.rou_w/self.rou_s -1
        self.C2=(self.area*(self.L-self.L_w))*(self.rou_w/self.rou_s)*self.const[4] +self.const[3]*self.area*self.L_w
        

    def dMw(self):

        Mw=self.W_sp-self.W_co+self.W_sg
        # W_sg needs to be introduced 
        return Mw
    
    def dLw(self):

        Lw=(1/(self.area*self.L_w))*(self.area*(self.L-self.L_w)*self.Kp5-self.C2/self.C1)*self.dPr()+self.W_sp/self.C1+self.W_co/self.C2\
        self.W_co=(self.C2*self.dPr()-self.W_sr-self.W_sp)/self.C1

        return Lw

    def dPr(self):

        """ constants
         Xp1,Xp2,Xp3,Xp4,Xp5,Xp6,Xp7
           need to be initialized in the constructor
        """

        self.meu_f=  self.Xp1+self.const[0]*self.P_pr
        self.meu_g=  self.Xp2+self.const[1]*self.P_pr
        self.meu_fg= self.Xp3+self.const[2]*self.P_pr
        self.hf=     self.Xp6+self.const[5]*self.P_pr
        self.hfg=    self.Xp7+self.const[6]*self.P_pr


        numerator=self.Q+self.W_sr*(self.P_pr*self.meu_g/(self.J*self.C1)+self.hfg/self.C1)\
        +self.W_sp*(self.h_sp-self.hf+self.hfg/self.C1+self.P_pr*(self.meu_g+self.meu_fg)/(self.J*self.C1))

        denominator=self.M_w*(self.const[5]+self.const[1]*self.P_pr/self.J)+self.M_s*self.const[1]*self.P_pr/self.J\
        -self.V_w/self.J +self.C2/self.C1*(self.hfg+self.P_pr*(self.meu_f+self.meu_fg)/self.J)

        Pr=numerator/denominator

        return Pr
