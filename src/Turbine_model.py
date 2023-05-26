import numpy as np
import matplotlib.pyplot as plt 
import threading
from pandas import ExcelWriter
import scipy as sp
from Condenser_Model import Condenser
from UTSG_MODEL import UTSG

"""
                         --------" ASTARA --a Nuclear Power Plant simulator "--------                                             
                                           
                        -------------------PRESSURIZER MODEL ------------------------

PROGRAMMER:

EBNY WALID AHAMMED 
Undergrad Student (Level 4 term 1)
Dept of Nuclear Engineering 
University of Dhaka                                          



"""
class Turbine():

    def __init__(self,Throttle_valve:object,moisture_seperator:object
                 ,Reheater:object):
        self.Throttle_valve=Throttle_valve
        self.moisture_seperator=moisture_seperator
        self.Reheater=Reheater
    
class Throttle_valve():
    def __init__(self,Area_main:float,Area_secondary:float,
                co_efficient_main:float,co_efficient_secondary:float,
                pos_main_valve:float,pos_second_valve:float,enthalpy_at_main_throttle_valve:float):
        
        """ constants """
        self.A_m=Area_main
        self.A_s=Area_secondary
        self.Cf_m=co_efficient_main
        self.Cf_s=co_efficient_secondary
        self.W_utsg=UTSG.Wst
        self.h_s=enthalpy_at_main_throttle_valve

        """ control elements """
        if ((pos_main_valve >100 and pos_main_valve<0) or( pos_second_valve>100 and pos_second_valve<0)):
            raise ValueError ("opening percentage of the  Main and secondarythrottle valve must stay in between 0 and 100")
        else:
            self.Pos_main_v=pos_main_valve
            self.pos_second_v=pos_second_valve

        # position of the main and second valve means how much the percentage of the valve is opened 
        #Here is Pos_main_V=100 then it's fully open
    def _Wmain(self):
        self.W_m=self.Cf_m*self.A_m*self.W_utsg*self.Pos_main_v
        self.W_2nd=self.Cf_s*self.A_s*self.W_utsg*self.pos_second_v

class Nozzle_chest(Throttle_valve):

    def __init__(self,Effective_volume_of_nozzle_Chest:float,steam_pressure_chest:float,
                 steam_density_chest:float,nozzle_chest_enthalpy:float,
                 Kc_hp:float,Callender_const1:float,Callender_const2:float):
        
        """ k1 and k2 are constants in the Callender’s emperical
            equation relating pressure, density and enthalpy of superheated steam
            """
        
        self.W_main=Throttle_valve.W_m
        self.h_sd=Throttle_valve.h_s
        self.rou_c=steam_density_chest
        self.k1=Callender_const1
        self.k2=Callender_const2
        self.Kc_hp=Kc_hp

        #variables 
        self.P_c=steam_pressure_chest
        self.V_c=Effective_volume_of_nozzle_Chest
        self.h_c=nozzle_chest_enthalpy
        
        
    def _Whp(self):
        self.Whp1=self.Kc_hp*np.sqrt(self.P_c*self.rou_c-Reheater.P*HighPressureTurbine.rou_exit)
 
    def _dh_c(self):
        Dhc=(((self.W_main*self.h_sd-self.Whp1*self.h_c)/(self.rou_c*self.V_c))+
             (self.P_c/self.rou_c**2)*self._drou_c())/(1-self.k1)
        
        return Dhc

    def _drou_c(self):
        Drou_c=(self.W_main-self.Whp1)/self.V_c
        return Drou_c
    
    def _P_c(self):
        self.P_c=self.rou_c*(self.k1*self.h_c-self.k2)


class Moisture_seperator():
    def __init__(self):
        pass


class Reheater():
    def __init__(self,pressure:float):
        self.P=pressure
        pass

class HighPressureTurbine():
    def __init__(self,exit_steam_density:float) :
        self.rou_exit=exit_steam_density
        pass

class LowPressureTurbine():
    def __init__(self) :
        pass
