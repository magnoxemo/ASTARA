import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import threading

from Turbine_model import Turbine
from Reactor_MODEL import Reactor
from Condenser_Pump_model import Condenser_pump


"""
Project NAME:
                         ------ ASTARA --a Nuclear Power Plant simulator ------- 
PROGRAMMER:

EBNY WALID AHAMMED 
Undergrad Student (Level 4 term 1)
Dept of Nuclear Engineering 
University of Dhaka

"""

class Condenser():

    def __init__(self,water_mass_inside_of_the_condenser:float,water_droplet_rate_hot_well:float,
                 water_condensatin_flow_rate:float,outlet_flow_rate_LPF:float,time_constant_condensation_process:float,
                 steam_flow_rate_into_condenser:float,oulet_water_enthalpy:float,enthalpy_of_saturated_water:float,
                 enthalpy_of_the_flow_entering:float,latent_enthalpy_of_water:float,water_flow_rate_to_condenser:float):
        
        """initials conditions + variables """

        self.Mw=water_mass_inside_of_the_condenser
        self.W1=water_droplet_rate_hot_well               #from turbine
        self.W2=steam_flow_rate_into_condenser
        self.W3=water_condensatin_flow_rate
        self.Wo=outlet_flow_rate_LPF                     #f toeed water heater 
        self.Wi=water_flow_rate_to_condenser
        self.ho=oulet_water_enthalpy

        """constants """
        self.tau_co=time_constant_condensation_process
        self.hfg=latent_enthalpy_of_water
        self.hi=enthalpy_of_the_flow_entering
        self.hf=enthalpy_of_saturated_water


    def dMw(self):

        Mw=self.W1-self.Wo+self.W3

        return Mw
    
    def dW3(self):

        W3=(self.W2-W3)/self.tau_co

        return W3
    
    def  dho(self):

        ho=(self.W1+self.W2)*(self.hf-self.ho)/self.Mw

        return ho
    
    def consitutive_eq(self):

        self.W1=self.Wi*(1-(self.hi-self.hf)/self.hfg)
        self.W2=(self.hi-self.hf)/self.hfg
