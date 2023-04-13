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
                secondaryFlow_area_in_the_U_tube_region:float,
                area_of_the_drum_water_section:float,pressure_drop_coefficient_in_the_recirculating_loop:float,
                Steam_valve_coefficient:float,Specific_heat_capacity_of_the_metaltubes:float,
                Specific_heat_capacity_of_primary_subcooled_region_fluid:float,
                Avg_enthalpy_boiling_region:float,Saturated_enthal_water:float,latent_enthal_water:float,
                Ex_enthal_boiling_region:float,water_level_inUTSG:float,
                Metal_mass_in_metal_node1:float,Metal_mass_in_metal_node2:float,
                Mass_water_in_the_the_primarynode1:float,Mass_water_in_the_the_primarynode2:float,
                Mass_water_in_the_the_primarynode3:float,Mass_water_in_the_the_primarynode4:float,
                Mass_of_water_in_t_h_e_inlet_plenum:float,Steam_generator_pressure:float,
                Heat_transfer_area_from_the_Utubes_to_the_secondary_side_in_the_subcooled:list,
                Heat_transfer_area_from_the_Utubes_to_the_secondary_side_in_boiling_region:float,
                Heat_transfer_areas_from_the_primary_side_to_the_Utubes_in_node1:float,
                Heat_transfer_areas_from_the_primary_side_to_the_Utubes_in_node2:float,Downcomer_temperature:float,
                Drum_water_temperature:float,Metal_tube_temperature_node1:float,Metal_tube_temperature_node2:float,
                Metal_tube_temperature_node3:float,Metal_tube_temperature_node4:float,
                Primary_coolant_temperatures_node1:float,Primary_coolant_temperatures_node2:float,
                Primary_coolant_temperatures_node3:float,Primary_coolant_temperatures_node4:float,
                Coolant_temperature_in_the_inlet_plenum:float,Coolant_temperature_in_the_outlet_plenum:float,
                Saturated_temperature_of_the_water_and_steam_in_the_UTSG:float,
                Heat_transfer_coefficient_from_the_primary_to_metal:float,
                Heat_transfer_coefficient_from_the_metal_to__subcooled:float,
                Heat_transfer_coefficient_from_the_metal_to_boiling_regions:float,
                Volume_of_the_drum_section:float, Specific_volume_of_the_saturated_water:float,
                Specific_volume_of_the_saturated_steam:float,Volume_of_the_riser_region:float,
                Steam_flow_rate:float,Constant_parameters_of_the_water_property_equations:list,
                steamquality_leaving_the_boiling_region:float,
                Average_density_of_the_fluid_in_the_boiling_region:float,
                Density_of_the_saturated_steam:float,Density_of_the_fluid_riser_region:float) -> None:           
        
                 #inlet temp of primary water= reactor's hot leg temp 
                 #inlet_temp_feed water=exit temp of the condenser

        

        pass
