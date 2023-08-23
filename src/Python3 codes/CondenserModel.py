import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

class Condenser():

    def __init__(self,water_mass_inside_of_the_condenser:float,water_droplet_rate_hot_well:float,
                 water_condensatin_flow_rate:float,outlet_flow_rate_LPF:float,time_constant_condensation_process:float,
                 steam_flow_rate_into_condenser:float,oulet_water_enthalpy:float,TemperatureIndise:float,
                 enthalpy_of_the_flow_entering:float,water_flow_rate_to_condenser:float):
        
        """ there are also some rooms to add more features for example 
            3rd loop water exit and inlet temperature 
            cooling tower dynamic modeling here. But for now the working modeling of
            the condenser is done  """
        
        """initials conditions + variables """

        self.Mw=water_mass_inside_of_the_condenser
        self.W1=water_droplet_rate_hot_well               
        self.W2=steam_flow_rate_into_condenser
        self.W3=water_condensatin_flow_rate
        self.Wo=outlet_flow_rate_LPF                      #feed water heater 
        self.Wi=water_flow_rate_to_condenser              #from turbine
        self.ho=oulet_water_enthalpy
        self.Temp=TemperatureIndise


        """constants """
        self.time_const=time_constant_condensation_process

        self.hf=PropsSI("H","T",self.Temp,"Q",0,"water")
        self.hg=PropsSI("H","T",self.Temp,"Q",1,"water")
        self.hfg=self.hg-self.hf

        self.hi=enthalpy_of_the_flow_entering


    def DMw(self):

        dtdMw=self.Wi*(1-(self.hi-self.hf)/self.hfg)-self.Wo+self.W3

        return dtdMw
    
    def DW3(self):

        dtdW3=(self.self.Wi*(self.hi-self.hf)/self.hfg-self.W3)/self.time_const

        return dtdW3
    
    def  Dho(self):

        dtdho=(self.Wi*(1-(self.hi-self.hf)/self.hfg)+self.W3)*(self.hf-self.ho)/self.Mw

        return dtdho
    
    def consitutive_eq(self):
    
        '''this equation is just to monitor the  
        
                1.condensation
                2.Water doplet from the vapor
                3.Vapor condensation rate when comes in contact with the 3rd loop water '''

        self.W1=self.Wi*(1-(self.hi-self.hf)/self.hfg)
        self.W2=self.Wi*(self.hi-self.hf)/self.hfg


    def integrator(self,function,argsforfunction:list,intitial_cond,time_step):
        l=len(argsforfunction)

        if l==0:
            return function()*time_step+intitial_cond
        elif l==1:
            arg1=argsforfunction[0]
            return function(arg1)*time_step+intitial_cond  
        elif l==2:
            arg1=argsforfunction[0]
            arg2=argsforfunction[1]
            return function(arg1,arg2)*time_step+intitial_cond
        elif l==3:
            arg1=argsforfunction[0]
            arg2=argsforfunction[1]
            arg3=argsforfunction[2]
            return function(arg1,arg2,arg3)*time_step+intitial_cond  
        else:
            raise   AttributeError("agrs in your differential function were not correct! Fix them")
        
