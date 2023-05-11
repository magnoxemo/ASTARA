import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sp 
from Reactor_MODEL import Reactor
from UTSG_MODEL import UTSG
from reactor_pump import pump

class Pressurizer():
    def __init__(self,area:float,water_enthalpy:float,latent_heat_of_vaporization:float,
                 water_enthalpy_spray_nozzle:float,heat_conversion_factor:float,
                 effective_length_of_pressurizer:float,water_level:float,
                 mass_steam:float,mass_water:float,pressure_at_pressurizer:float,
                 condensation_flow_rate:float,spray_flow_rate:float,in_and_out_flow_rate:float,
                 water_vol_in_pressurizer:float,specific_vol_water:float,specific_vol_steam:float,
                 steam_density:float,water_density:float):
        
        #design parameters
        self.area=area
        self.L=effective_length_of_pressurizer
        self.J=heat_conversion_factor

        #constants 
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
        self.P_w=water_vol_in_pressurizer
        self.P_pr=pressure_at_pressurizer
        self.W_co=condensation_flow_rate

    #math functions and model will be updated here 
