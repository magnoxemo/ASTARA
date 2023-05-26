import numpy as np
import matplotlib.pyplot as plt 
import threading
from pandas import ExcelWriter
import scipy as sp
from Condenser_Model import Condenser
from UTSG_MODEL import UTSG


class Turbine():

    def __init__(self,Throttle_valve:object,moisture_seperator:object
                 ) -> None:
        Throttle_valve=Throttle_valve
        moisture_seperator=moisture_seperator
    
class Throttle_valve():
    def __init__(self,Area_main:float,Area_secondary:float,
                co_efficient_main:float,co_efficient_secondary:float,
                pos_main_valve:float,pos_second_valve:float):
        
        """ constants """
        self.A_m=Area_main
        self.A_s=Area_secondary
        self.Cf_m=co_efficient_main
        self.Cf_s=co_efficient_secondary
        self.W_utsg=UTSG.Wst

        """ control elements """
        if ((pos_main_valve >100 and pos_main_valve<0) or( pos_second_valve>100 and pos_second_valve<0)):
            raise ValueError ("opening percentage of the  Main and secondarythrottle valve must stay in between 0 and 100")
        else:
            self.Pos_main_v=pos_main_valve
            self.pos_second_v=pos_second_valve

        # position of the main and second valve means how much the percentage of the valve is opened 
        #Here is Pos_main_V=100 then it's fully open
    def Wmain(self):
        self.W_main=self.Cf_m*self.A_m*self.W_utsg*self.Pos_main_v
        self.W_2nd=self.Cf_s*self.A_s*self.W_utsg*self.pos_second_v
        
class Moisture_seperator():
    def __init__(self) -> None:
        pass
