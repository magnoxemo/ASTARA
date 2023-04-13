import numpy as np
import math
import matplotlib.pyplot as plt 

"""
Project NAME:
                         --------" ASTARA --a Nuclear Power Plant simulator "-------- 
PROGRAMMER:

EBNY WALID AHAMMED 
Undergrad Student (Level 4 term 1)
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
      will be two different objects but there will be a function for calculations and data interchange 

from model prespective:
    1.The whole reactor will be considered as a uniform fuel rod 
    2.The core will be divided in different parts (15+2 lumps )
    3.A graph of the temperature distribution will be shown 
    4.Modeling will follow the same assumption as the desetation 
    5.for solving the differential equation Runge Kutta 4th order method will be used ()

     
"""

class Reactor ():
    P0=1

    def __init__(self,Precursor_concentration:float,Area:float,
                 Coolant_heat_capacity:float,Fuel_heat_capacity:float,
                 Fission_power_fraction:float,Average_heat_transfer_factor:float
                 ,Coolant_mas_flow_rate:float,Coolant_mass_in_twofluid_nodes:float,
                 Clod_leg_water_mass:float,fuel_nodal_mass:float,
                 Hot_leg_water_mass:float,Lower_plenum_water_mass:float,
                 Coolant_node_mass:float,upper_plenum_water_mass:float,
                 Reactor_core_power:float,Cold_leg_temp:float,
                 hot_leg_tem:float,Fluid_temp_in_lower_plenum:float,
                 Fluid_temp_in_upper_plenum:float,Fuel_temp_in_each_node:list,
                 Moderatro_temp_in_nodes:list,Steam_gen_inlet_temp:float,
                 Reactivity_coefficient_coolant:float,Reactivity_coefficient_fuel:float,beta:list,
                 neutron_group_const:list,neutron_generation_time:float,
                 Total_reactivity:float,external_reactivity:float,
                 TempF0_initial_cond:list,TempM0_initial_cond:list) -> None:
        
        """local variables with appropriate notations"""
        self.C=Precursor_concentration
        self.A=Area
        self.h=Average_heat_transfer_factor
        self.C_pc=Coolant_heat_capacity
        self.C_pf=Fuel_heat_capacity
        self.F_r=Fission_power_fraction
        self.M=Coolant_mas_flow_rate
        self.M_c=Coolant_mass_in_twofluid_nodes
        self.M_cl=Clod_leg_water_mass
        self.M_f=fuel_nodal_mass
        self.M_hl=Hot_leg_water_mass
        self.M_lp=Lower_plenum_water_mass
        self.M_mo=Coolant_node_mass
        self.M_up=upper_plenum_water_mass
        self.P=Reactor_core_power
        self.T_cl=Cold_leg_temp
        self.T_fuel_node1=Fuel_temp_in_each_node[0]
        self.T_fuel_node2=Fuel_temp_in_each_node[1]
        self.T_fuel_node3=Fuel_temp_in_each_node[2]
        self.T_hl=hot_leg_tem
        self.T_lp=Fluid_temp_in_lower_plenum
        self.T_up=Fluid_temp_in_upper_plenum
        self.T_moderator_node1=Moderatro_temp_in_nodes[0]
        self.T_moderator_node2=Moderatro_temp_in_nodes[1]
        self.T_moderator_node3=Moderatro_temp_in_nodes[2]
        self.T_moderator_node4=Moderatro_temp_in_nodes[3]
        self.T_moderator_node5=Moderatro_temp_in_nodes[4]
        self.T_moderator_node6=Moderatro_temp_in_nodes[5]
        self.T_po=Steam_gen_inlet_temp
        self.alpha_c=Reactivity_coefficient_coolant
        self.alpha_f=Reactivity_coefficient_fuel
        self.beta=np.array(beta)
        self.Beta_t=np.sum(beta)
        self.lemda=np.array(neutron_group_const)
        self.T_neutron_g=len(self.lemda)
        self.NGT=neutron_generation_time
        self.Ro=Total_reactivity
        self.Ro_ex=external_reactivity
        self.Tf0=np.array(TempF0_initial_cond)
        self.TM0=np.array(TempM0_initial_cond)

        

    def constitutive_cal(self):
        sum_mod1=self.T_moderator_node1+self.T_moderator_node2+self.T_moderator_node3
        +self.T_moderator_node4+self.T_moderator_node5+self.T_moderator_node6

        sum_f1=self.T_fuel_node1+self.T_fuel_node2+self.T_fuel_node3

        sum_mod0=np.sum(self.TM0)
        sum_f0=np.sum(self.Tf0)


        self.Tau_cl=self.M_cl/self.M
        self.Tau_hl=self.M_hl/self.M
        self.Tau_C=self.M_c/(2*self.M)
        self.Tau_lp=self.M_lp/self.M
        self.Tau_up=self.M_up/self.M

        self.LAMDA =self.Beta_t/np.sum(self.beta/self.lemda)
        self.Ro=self.Ro_ex+self.alpha_c*(sum_mod1-sum_mod0)/6+self.alpha_f*(sum_f1-sum_f0)/3

    def dT_lp(self):

        Tlp=(self.T_cl-self.T_lp)/self.Tau_lp

        return Tlp
    
    def dT_up(self):

        Tup=(self.T_moderator_node6-self.T_up)/self.Tau_up

        return Tup
    
    def dT_hl(self):

        Thl=(self.T_up-self.T_hl)/self.Tau_hl

        return Thl
    

    def dT_cl(self): #---connect to feed water system

        Tcl=(self.T_po-self.T_cl)/self.Tau_cl

        return Tcl
    

    def dTf1(self):

        Tf1=self.F_r*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(self.T_moderator_node1-self.T_fuel_node1)/(self.M*self.C_pf)

        return Tf1
    
    def dTf2(self):

        Tf2=self.F_r*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(self.T_moderator_node3-self.T_fuel_node2)/(self.M*self.C_pf)

        return Tf2
    
    def dTf3(self):

        Tf3=self.F_r*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(self.T_moderator_node5-self.T_fuel_node3)/(self.M*self.C_pf)

        return Tf3
    
    def dTmol1(self):

        TMol1= (1-self.F_r)*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(-self.T_moderator_node1+self.T_fuel_node1)/(self.M*self.C_pf)
        +(self.T_lp-self.T_moderator_node1)/self.Tau_C

        return TMol1
        pass 
    def dTmol2(self):

        TMol2= (1-self.F_r)*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(-self.T_moderator_node1+self.T_fuel_node1)/(self.M*self.C_pf)
        +(self.T_moderator_node1-self.T_moderator_node2)/self.Tau_C
        
        return TMol2
    
    def dTmol3(self):

        TMol3= (1-self.F_r)*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(-self.T_moderator_node3+self.T_fuel_node2)/(self.M*self.C_pf)
        +(self.T_moderator_node2-self.T_moderator_node3)/self.Tau_C

        return TMol3
        pass 
    def dTmol4(self):

        TMol4= (1-self.F_r)*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(-self.T_moderator_node3+self.T_fuel_node2)/(self.M*self.C_pf)
        +(self.T_moderator_node3-self.T_moderator_node4)/self.Tau_C
        
        return TMol4
    
    def dTmol5(self):

        TMol5= (1-self.F_r)*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(-self.T_moderator_node5+self.T_fuel_node3)/(self.M*self.C_pf)
        +(self.T_moderator_node4-self.T_moderator_node5)/self.Tau_C

        return TMol5
 
    def dTmol6(self):

        TMol6= (1-self.F_r)*self.P0*(self.P/self.P0)/(self.M*self.C_pf)
        +self.h*self.A*(-self.T_moderator_node5+self.T_fuel_node3)/(self.M*self.C_pf)
        +(self.T_moderator_node5-self.T_moderator_node6)/self.Tau_C
        
        return TMol6
    
    def step_integator(self,function,condition:list,stepsize):

        dt=stepsize
        y0=condition[0]

        return y0+function()*dt
    
"""reactor modeling
things remains:
1.Have to couple this with the steam generator 
2.Have to couple this with the feed water pump 
3.Have to couple this with the reactor controller 
4.Have to write the function-- all the calculation initiator 
"""
