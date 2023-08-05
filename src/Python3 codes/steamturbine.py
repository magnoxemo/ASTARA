import numpy as np
import matplotlib.pyplot as plt 
import threading 
import scipy as sp

class inlet_plenum():
    
    def __init__(self,Temperature:float):
        """
        
        Mass     ---> stays in the pipe under steady state conditions
        Flow_rate---> hot_leg_flow_rate
        Theta    ---> temperature of the hot leg water 
        
        """
        self.Mass=10000		 
        self.Flow_rate=4964.96
        self.density=732.134 
        self.Theta=400
               
        self.Temperature=Temperature
        self.time_const=self.Flow_rate/self.Mass
        
    def DT_pi(self,reactor:object):


        dtpi=(self.Theta-self.Temperature)*self.time_const
        
        """
        self.theta=reactor.T_hotleg
        self.theta should be coupled with the reactor 
        it will be reactor.T_hotleg
        
        """

        return dtpi
    
    def integrator(self,function,intitial_cond,time_step):

        return function()*time_step+intitial_cond



class primary_lump():

    def __init__(self,PrimaryLumpTemperature:list,MetalLumpTemperature:list,outletTemperature:float):

        self.number_of_utube=3383
        self.inner_diameter=0.019685
        self.outer_diameter=0.022225
        self.density=732.134 
        self.length=10.831712
        self.frist_lump_length=1.05017116
        self.heat_capacity_1=5819.652

        self.Primary_side_flim_conductance=25563
        self.Tube_metal_conductance=12263.68
        self.Tube_metal_conductance_subcool=11186.216
        self.Tube_metal_conductance_boiling=34068

        self.Wpi=4964.96   #hot_leg_flow_rate
        self.wfi=470.226   #Trubine_outlet
        self.W1=2349.45    #SFSL
        self.W2=self.W1
        self.W3=self.W1
        
        self.second_lump_length=self.length-self.frist_lump_length
        self.Ap=np.pi*self.inner_diameter**2*self.number_of_utube/4
        self.Mp=self.density*self.Ap*self.length
        self.Mp1=self.Mp*(self.frist_lump_length/self.length)
        self.Mp2=self.Mp*(self.second_lump_length/self.length)
        self.Mp3=self.Mp2
        self.Mp4=self.Mp1

        self.Sm=np.pi*self.length*self.outer_diameter*self.number_of_utube
        self.Sm1=self.Sm*self.frist_lump_length/self.length
        self.Sm2=self.Sm*self.second_lump_length/self.length
        self.Sm3=self.Sm2
        self.Sm4=self.Sm1

        self.Spm1=self.Sm1*self.inner_diameter/self.outer_diameter
        self.Spm2=self.Sm2*self.inner_diameter/self.outer_diameter
        self.Spm3=self.Spm2
        self.Spm4=self.Spm1

        self.Pr1=self.Spm1/self.frist_lump_length
        self.Pr2=self.Spm2/self.second_lump_length

        """
        initial conditions (as this constructor will only run once)
        PrimaryLumpTemperature will carry the initial condition for 4 primary lumps
        MetalLumpTemperature will carry the initial condition for 4 metal lumps as per model
        MetalLumpTemperature will be initialized in the Metal lump object
        """ 
        if len(PrimaryLumpTemperature)!=4:
            raise ValueError(" Initial condition error!")
        else:
            self.Tp1=PrimaryLumpTemperature[0]
            self.Tp2=PrimaryLumpTemperature[1]
            self.Tp3=PrimaryLumpTemperature[2]
            self.Tp4=PrimaryLumpTemperature[3]

        if len(MetalLumpTemperature)!=4:
            raise ValueError(" Initial condition error!")
        else:
            self.Tm1=MetalLumpTemperature[0]
            self.Tm2=MetalLumpTemperature[1]
            self.Tm3=MetalLumpTemperature[2]
            self.Tm4=MetalLumpTemperature[3]

        self.Tpo=outletTemperature
         
    
    def DTp1(self,inlet_plenum:object):
        """inlet_plenum"""
        dtdTp1=self.Wpi*(inlet_plenum.Temperature-self.Tp1)/(self.density*self.Ap*self.frist_lump_length)\
              +self.Primary_side_flim_conductance*self.Spm1*(self.Tm1-self.Tp1)\
              /(self.Mp1*self.heat_capacity_1)
        return dtdTp1
    
    def DTp2(self,sub_cool_region:object):
        '''DLs1 method will be called from the sub_cool_region class '''
        dtdTp2=self.Wpi*(self.Tp1-self.Tp2)/(self.density*self.Ap*self.second_lump_length)\
              +self.Primary_side_flim_conductance*self.Spm2*(self.Tm2-self.Tp2)\
              /(self.Mp1*self.heat_capacity_1)+(self.Tp2-self.Tp1)*sub_cool_region.DLs1()/\
              self.second_lump_length
        return dtdTp2
    
    def DTp3(self):
        dtdTp3=self.Wpi*(self.Tp2-self.Tp3)/(self.density*self.Ap*self.frist_lump_length)\
              +self.Primary_side_flim_conductance*self.Spm2*(self.Tm3-self.Tp3)\
              /(self.Mp1*self.heat_capacity_1)       
        return dtdTp3
    
    def DTp4(self,sub_cool_region:object):
        '''DLs1 method will be called from the sub_cool_region class '''
        dtdTp4=self.Wpi*(self.Tp1-self.Tp2)/(self.density*self.Ap*self.second_lump_length)\
              +self.Primary_side_flim_conductance*self.Spm1*(self.Tm4-self.Tp4)\
              /(self.Mp1*self.heat_capacity_1)+(self.Tp3-self.Tp4)*sub_cool_region.DLs1()/\
              self.frist_lump_length
        return dtdTp4   
     
    def integrator(self,function,intitial_cond,time_step):
        return function()*time_step+intitial_cond
