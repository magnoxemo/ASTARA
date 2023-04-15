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
                Tube_D_inside:float,Tube_height:float,UTSG_height:float,downcomer_length:float,
                secondaryFlow_area_in_the_U_tube_region:float,
                area_of_the_drum_water_section:float,pressure_drop_coefficient_in_the_recirculating_loop:float,
                Steam_valve_coefficient:float,Specific_heat_capacity_of_the_metaltubes:float,
                Specific_heat_capacity_of_primary_fluid:float,Specific_heat_capacity_of_subcooled_region_fluid:float,
                Avg_enthalpy_boiling_region:float,Saturated_enthal_water:float,latent_enthal_water:float,
                Ex_enthal_boiling_region:float,constant:list,water_level_inUTSG:float,subcooled_length:float,
                Metal_mass_in_metal_node1:float,Metal_mass_in_metal_node2:float,
                Mass_water_in_the_the_primarynode1:float,Mass_water_in_the_the_primarynode2:float,
                Mass_water_in_the_the_primarynode3:float,Mass_water_in_the_the_primarynode4:float,
                Mass_of_water_in_t_h_e_inlet_plenum:float,Steam_generator_pressure:float,
                Heat_transfer_area_from_the_Utubes_to_the_secondary_side_in_the_subcooled:float,
                Heat_transfer_area_from_the_Utubes_to_the_secondary_side_in_boiling_region:float,
                Heat_transfer_areas_from_the_primary_side_to_the_Utubes_in_node1:float,
                Heat_transfer_areas_from_the_primary_side_to_the_Utubes_in_node2:float,Downcomer_temperature:float,
                Drum_water_temperature:float,Metal_tube_temperature_node1:float,Metal_tube_temperature_node2:float,
                Metal_tube_temperature_node3:float,Metal_tube_temperature_node4:float,
                Primary_coolant_temperatures_node1:float,Primary_coolant_temperatures_node2:float,
                Primary_coolant_temperatures_node3:float,Primary_coolant_temperatures_node4:float,
                Coolant_temperature_in_the_inlet_plenum:float,Coolant_temperature_in_the_outlet_plenum:float,
                Saturated_temperature_of_the_water_and_steam_in_the_UTSG:float,Heat_transfer_coefficient_from_the_primary_to_metal:float,
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

        self.N=Number_of_Utube
        self.A_fs=secondaryFlow_area_in_the_U_tube_region
        self.A_dw=area_of_the_drum_water_section
        self.C1=pressure_drop_coefficient_in_the_recirculating_loop
        self.Cl=Steam_valve_coefficient
        self.Cm=Specific_heat_capacity_of_the_metaltubes
        self.Cp1=Specific_heat_capacity_of_primary_fluid
        self.Cp2=Specific_heat_capacity_of_subcooled_region_fluid
        self.H_b=Avg_enthalpy_boiling_region
        self.H_f=Saturated_enthal_water
        self.H_fg=latent_enthal_water
        self.H_ex=Ex_enthal_boiling_region

        self.K1=constant[0]
        self.K2=constant[1]
        self.K3=constant[2]
        self.K4=constant[3]
        self.K5=constant[4]
        self.K6=constant[5]

        self.L=UTSG_height
        self.L_d=downcomer_length
        self.L_dw=water_level_inUTSG
        self.L_s1=subcooled_length
        self.M_m1=Metal_mass_in_metal_node1
        self.M_m2=Metal_mass_in_metal_node2
        self.M_p1=Mass_water_in_the_the_primarynode1
        self.M_p2=Mass_water_in_the_the_primarynode2
        self.M_p3=Mass_water_in_the_the_primarynode3
        self.M_p4=Mass_water_in_the_the_primarynode4
        self.M_pi=Mass_of_water_in_t_h_e_inlet_plenum
        self.P=Steam_generator_pressure
        self.P_r1=Tube_D_inside
        self.P_r2=Tube_D_outside
        self.S_ms2=Heat_transfer_area_from_the_Utubes_to_the_secondary_side_in_boiling_region
        self.S_ms1=Heat_transfer_area_from_the_Utubes_to_the_secondary_side_in_the_subcooled
        self.S_pm1=Heat_transfer_areas_from_the_primary_side_to_the_Utubes_in_node1
        self.S_pm2=Heat_transfer_areas_from_the_primary_side_to_the_Utubes_in_node2

        self.T_d=Downcomer_temperature
        self.T_dw=Drum_water_temperature
        self.T_m1=Metal_tube_temperature_node1
        self.T_m2=Metal_tube_temperature_node2
        self.T_m3=Metal_tube_temperature_node3
        self.T_m4=Metal_tube_temperature_node4
        self.T_p1=Primary_coolant_temperatures_node1
        self.T_p2=Primary_coolant_temperatures_node2
        self.T_p3=Primary_coolant_temperatures_node3
        self.T_p4=Primary_coolant_temperatures_node4
        self.T_pi=Coolant_temperature_in_the_inlet_plenum
        self.T_po=Coolant_temperature_in_the_outlet_plenum
        self.T_sat=Saturated_temperature_of_the_water_and_steam_in_the_UTSG
        self.U_pm=Heat_transfer_coefficient_from_the_primary_to_metal
        self.U_ms1=Heat_transfer_coefficient_from_the_metal_to__subcooled
        self.U_ms2=Heat_transfer_coefficient_from_the_metal_to_boiling_regions
        self.V_dr=Volume_of_the_drum_section
        self.V_g=Specific_volume_of_the_saturated_steam
        self.V_f=Specific_volume_of_the_saturated_water
        self.V_fg=self.V_g-self.V_f
        self.V_r=Volume_of_the_riser_region
        self.Wst=Steam_flow_rate
        self.X1=Constant_parameters_of_the_water_property_equations[0]
        self.X2=Constant_parameters_of_the_water_property_equations[1]
        self.X3=Constant_parameters_of_the_water_property_equations[2]
        self.X4=Constant_parameters_of_the_water_property_equations[3]
        self.X5=Constant_parameters_of_the_water_property_equations[4]
        self.X6=Constant_parameters_of_the_water_property_equations[5]
        self.Xe=steamquality_leaving_the_boiling_region
        self.Rou_b=Average_density_of_the_fluid_in_the_boiling_region
        self.Rou_g=Density_of_the_saturated_steam
        self.Rou_r=Density_of_the_fluid_riser_region


    def constitutive_eq(self):


        self.H_f=self.X3+self.K3*self.P
        self.H_fg=self.x4+self.K4*self.P
        self.H_b=self.H_f+self.Xe*self.H_fg/2
        self.H_ex=self.H_f+self.Xe*self.H_fg
        self.L_s2=self.L-self.L_s1

        self.T_sat=self.X5+self.K5*self.P
        self.V_f=self.X1+self.K1*self.P
        self.V_fg=self.X2+self.K2*self.P
        self.Wst=self.C1*self.P
        self.Rou_b=1/(self.V_f+self.Xe*self.V_fg/2)
        self.Rou_r=1/(self.V_f+self.Xe*self.V_fg)
        self.Rou_g=self.X6+self.K6*self.P

        """every variable I calculate here don't need to be on the constructor
           and    will update it later """
    

