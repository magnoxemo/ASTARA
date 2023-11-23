import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import numpy as np
from CoolProp.CoolProp import PropsSI


class Condenser():
    def __init__(self,steampressure:float,airpressure:float):

        self.volume=3
        self.UA=356.972e3
        self.Cp=4.2e3
        self.Ps=steampressure
        self.Pa=airpressure
        self.Rs=0.4615e3
        self.Ra=.287e3

        """ steam mass balance variables"""
        self.W_turbine=4
        self.W_otherthanturbine=10
        self.W_condensate=4
        self.W_steamairout=0

        self.T_steamin=600

        
        """ air zone """
        self.W_vaccumbreakvalve=0
        self.W_air=0
        self.W_steamgas=0
        self.W_draincondenser=0

        """ hot well water"""
        self.Hot_wellarea=0.2
        self.W_hotwell=110
        self.W_bubblingoxygen=10

        """ cold water """
        self.W_coldwater=107.881
        self.T_coldwater=273+60
        self.W_hotwater=1010
        self.T_hotwater=273+80


        self.Hs=PropsSI('H','T',self.T_steamin,'P',self.Ps,'water')
        self.Hcw=PropsSI('H','T',self.T_coldwater,'P',101325,'water')


    def DW_steam(self):

        self.R=(self.Pa*self.Rs)/(self.Pa*self.Ra+self.Ps*self.Rs)
        self.Wss=self.W_air*(1-self.R)
        dtdWsteam=self.W_turbine+self.W_otherthanturbine-self.W_condensate-self.Wss

        return dtdWsteam
    
    def DPs(self):
        dtdPs=self.Rs*(self.DW_steam()*self.T_steamin)/self.volume
        return dtdPs
    def Denthalpy(self):
        pass 

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
        
