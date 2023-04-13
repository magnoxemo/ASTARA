import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sp 

"""
Project NAME:
                         -------" ASTARA --a Nuclear Power Plant simulator "-------- 
PROGRAMMER:

EBNY WALID AHAMMED 
Undergrad Student (Level 4 term 1)
Dept of Nuclear Engineering 
University of Dhaka



This will be a  generalized model of PWR.This code is based on this desetation: 

    Lumped Parameter, State Variable Dynamic Models for U-tube Recirculation Type Nuclear Steam Generators
    Mohamed Rabie Ahmed Ali
    University of Tennessee - Knoxville

    ------------------------ UTSG ---------------------------
    Model type D

     
"""

class UTSG():

    def __init__(self,Number_of_Utube:int,Tube_D_outside:float,
                 Tube_D_inside:float,Tube_height:float,UTSG_height:float,
                 flow_area_in_tube:float,flow_area_in_downcomer:float,
                 flow_area_in_riser:float,flow_area_in_drum:float,riser_height:float,
                 Primary_water_flow_rate:float,Primary_water_volume:float,
                 specific_heat_primary_water:float,inlet_temp_pwater:float,
                 outlet_temp_pwater:float,P_avg_Pside:float,avg_density:float,
                 outlet_steam_flow_rate:float,steam_pressure:float,
                 steam_temp_at_Spressure:float,inlet_temp_feed_water:float,
                 avg_density_secondary_subcooled_water:float,E_heat_transfer_area:float,
                 FHTC_primary_water:float,FHTC_secondary_subcooled_water:float,FHTC_secondary_boiling_water:float,
                 Metal_Tube_conductivity:float) -> None:           
        
                 #inlet temp of primary water= reactor's hot leg temp 
                 #inlet_temp_feed water=exit temp of the condenser

        pass
