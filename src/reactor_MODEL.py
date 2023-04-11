import numpy as np
import math
import matplotlib.pyplot as plt 

"""
Project NAME:
                         --------" ASTARA "------a Nuclear Power Plant simulator 
PROGRAMMER:

EBNY WALID AHAMMED 
Undergrad Student (term 1 level 4)
Dept of Nuclear Engineering 
University of Dhaka



This will be a  generalized model of PWR.This code is based on this desetation: 
https://trace.tennessee.edu/cgi/viewcontent.cgi?article=4038&context=utk_gradthes 

Concepts and features: 

from coding prespective:
    1.Reactor will be an OBJECT so by adjusting the parameters any reactor can be modeled
    2.Multithreading will be done to solve the nonlinear differential eq
    3.Reactor constants will be initialized via the OBJECT's constructor 
    4.Reactor and Reactor controller ---(in that case control rod driving system)
      will be two different objects but there will be a friend function for calculations and data interchange 
from model prespective:
    1.The whole reactor will be considered as a uniform fuel rod 
    2.The core will be divided in different parts (15+2 lumps )
    3.A graph of the temperature distribution will be shown 
    4.Modeling will follow the same assumption as the desetation 
    5.for solving the differential equation Runge Kutta 4th order method will be used 
     
"""

class Reactor ():

    
    def __init__(self,Precursor_concentration:float,Area:float,
                 Coolant_heat_capacity:float,Fuel_heat_capacity:float,
                 Fission_power_fraction:float,Average_heat_transfer_factor:float
                 ,Coolant_mas_flow_rate:float,Coolant_mass_in_twofluid_nodes:float,
                 Clod_leg_water_mass:float,Number_of_nodes:int,total_fuel_mass:float,
                 Hot_leg_water_mass:float,Lower_plenum_water_mass:float,
                 Coolant_node_mass:float,upper_plenum_water_mass:float,
                 Reactor_core_power:float,Cold_leg_temp:float,
                 hot_leg_tem:float,Fluid_temp_in_lower_plenum:float,
                 Fluid_temp_in_upper_plenum:float,Fuel_temp_in_each_node:list,
                 Moderatro_temp_in_nodes:list,Steam_gen_inlet_temp:float,
                 Reactivity_coefficient_coolant:float,Reactivity_coefficient_fuel:float,
                 Total_delayed_neutron_group_fraction:float,
                 neutron_group_const:list,neutron_generation_time:float,
                 Total_reactivity:float,external_reactivity:float,timeconst_coldLeg:float,
                 timeconst_hotLeg:float,timeconst_moderatorNode:float,
                 timeconst_lower_plenum:float,timeconst_upper_plenum:float) -> None:
        
        """local variables with appropriate notations"""
        self.C=Precursor_concentration
        self.A=Area
        self.Node_num=Number_of_nodes
        self.C_pc=Coolant_heat_capacity
        self.C_pf=Fuel_heat_capacity
        self.F_r=Fission_power_fraction
        self.M=Coolant_mas_flow_rate
        self.M_c=Coolant_mass_in_twofluid_nodes
        self.M_cl=Clod_leg_water_mass
        self.M_f=total_fuel_mass/self.Node_num
        self.M_hl=Hot_leg_water_mass
        self.M_lp=Lower_plenum_water_mass
        self.M_mo=Coolant_node_mass
        self.M_up=upper_plenum_water_mass
        self.P=Reactor_core_power
        self.T_cl=Cold_leg_temp
        self.T_fuel_node=np.array(Fuel_temp_in_each_node)
        self.T_hl=hot_leg_tem
        self.T_lp=Fluid_temp_in_lower_plenum
        self.T_up=Fluid_temp_in_upper_plenum
        self.T_moderator_node=np.array(Moderatro_temp_in_nodes)
        self.T_po=Steam_gen_inlet_temp
        self.alpha_c=Reactivity_coefficient_coolant
        self.alpha_f=Reactivity_coefficient_fuel
        self.Beta_t=Total_delayed_neutron_group_fraction
        self.lemda=np.array(neutron_group_const)
        self.T_neutron_g=len(self.lemda)
        self.NGT=neutron_generation_time
        self.Ro=Total_reactivity
        self.Ro_ex=external_reactivity
        self.Tau_cl=timeconst_coldLeg
        self.Tau_hl=timeconst_hotLeg
        self.Tau_C=timeconst_moderatorNode
        self.Tau_lp=timeconst_lower_plenum
        self.Tau_up=timeconst_upper_plenum


        if len(self.T_fuel_node)==self.Node_num:

            pass

        else:

            raise AttributeError('Fuel node temperature missing')
        
        if len(self.T_moderator_node)==2*self.Node_num:
            pass
        else:
            raise AttributeError('Moderator node temperature missing')
