import numpy as np 
import matplotlib.pyplot as plt 
import scipy as sp
import threading as th

""" Reactor modeling """

class CoreInside():
    def __init__(self,TemperatureFuel:list,TemperatureModerator:list,Power:float,Precursor:float):
        self.diameter=36.2712
        self.core_height=43.8912
        self.CoolantDensity=732.18278
        self.area=5564.8921/3
        self.h=50537.48

        #group const
        self.Beta1=0.000209
        self.Beta2=0.001414
        self.Beta3=0.001414
        self.Beta4=0.00272 
        self.Beta5=0.00092
        self.Beta6=0.00689

        self.total_delayed_const=self.Beta1+self.Beta2+self.Beta3+self.Beta4\
            +self.Beta5+self.Beta6
        
        #decay const
        self.Lamda1=0.0125
        self.Lamda2=0.0308
        self.Lamda3=0.114
        self.Lamda4=0.307
        self.Lamda5=1.19
        self.Lamda6=3.19

        self.NGT=17.9*10**-6
        self.Alpha_m=-1.1192*10**-7
        self.Alpha_f=6.111*10**-6

        self.Fr=0.974       #fission power factor 
        self.Cpf=247.5052   #coolant conductivity 
        self.Cpc=5819.652   #fuel heat conductivity 
        #volumes 
        self.CoolantVcore=15.2911

        self.FuelMass=101032.711
        self.Wc=18899.68208

        """initial conditions"""
        self.NominalPower=Power 
        self.Precursor=Precursor
        self.PowerRatio=1
        self.ExternalReactivity=0  #control parameter will be used by control driving program 

        if len(TemperatureFuel)!=3:
            raise ValueError("Here should be three initial conditions!")
        else:
            self.Tf1=TemperatureFuel[0]
            self.Tf2=TemperatureFuel[1]
            self.Tf3=TemperatureFuel[2]
        if len(TemperatureModerator)!=6:
            raise ValueError('Here should be six initial conditions!')
        else:
            self.Tmo1=TemperatureModerator[0]
            self.Tmo2=TemperatureModerator[1]
            self.Tmo3=TemperatureModerator[2]
            self.Tmo4=TemperatureModerator[3]
            self.Tmo5=TemperatureModerator[4]
            self.Tmo6=TemperatureModerator[5]

        """ necessary calculations """
        """mass"""
        self.Mf=self.FuelMass/3 #as there will three lumps 
        self.Mmo=self.CoolantVcore*self.CoolantDensity /3 

        """time const"""
        self.time_constmo=self.Mmo/(self.Wc*2)

        self.Lamda=self.total_delayed_const/(self.Beta1/self.Lamda1+self.Beta2/self.Lamda2+self.Beta3/self.Lamda3+self.Beta4/self.Lamda4\
                                             +self.Beta5/self.Lamda5+self.Beta6/self.Lamda6)
        
    def Reacivity(self,FuelTemp:list,ModeratorTemp:list):
        
        self.Reactivity=self.ExternalReactivity+self.Alpha_f*(self.Tf1+self.Tf2+self.Tf3-sum(FuelTemp))/3\
                        +(self.Tmo1+self.Tmo2+self.Tmo3+self.Tmo4+self.Tmo5+self.Tmo6-sum(ModeratorTemp))
    
    def DPowerRatio(self):
        dtdPP0=(self.Reacivity()-self.total_delayed_const)*self.PowerRatio/self.NGT+self.Lamda*self.Precursor
        return dtdPP0
        
    def DPrecoursor(self):
        dtdc=self.total_delayed_const*self.PowerRatio/self.NGT-self.Lamda*self.Precursor
        return dtdc
        
    def DTf1(self):
        dtdTf1=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo1-self.Tf1)/(self.Mf*self.Cpf)
        return dtdTf1
        
    def DTf2(self):
        dtdTf2=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo3-self.Tf2)/(self.Mf*self.Cpf)
        return dtdTf2
        
    def DTf3(self):
        dtdTf3=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo1-self.Tf1)/(self.Mf*self.Cpf)
        return dtdTf3
    
    def DTm01(self,LowerPlenum:object):
        dtdTm01=(1-self.fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)+\
            self.h*self.area*(LowerPlenum.Tf1-self.Tmo1)+(self.Tlp-self.Tmo1)/self.time_constmo
        return dtdTm01 
        
    def DTm02(self):
        dtdTm02=(1-self.fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)+\
            self.h*self.area*(self.Tf1-self.Tmo1)+(self.Tmo1-self.Tmo2)/self.time_constmo
        return dtdTm02 

    def DTm03(self):
        dtdTm03=(1-self.fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)+\
            self.h*self.area*(self.Tf2-self.Tmo3)+(self.Tmo2-self.Tmo3)/self.time_constmo
        return dtdTm03
        
    def DTm04(self):
        dtdTm04=(1-self.fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)+\
            self.h*self.area*(self.Tf2-self.Tmo3)+(self.Tmo3-self.Tmo4)/self.time_constmo
        return dtdTm04
    
    def DTm05(self):
        dtdTm05=(1-self.fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)+\
            self.h*self.area*(self.Tf3-self.Tmo5)+(self.Tmo4-self.Tmo5)/self.time_constmo
        return dtdTm05
        
    def DTm06(self):
        dtdTm06=(1-self.fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)+\
            self.h*self.area*(self.Tf3-self.Tmo5)+(self.Tmo5-self.Tmo6)/self.time_constmo
        return dtdTm06
    
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
    
class Hotleg():
    def __init__(self,TempHotLeg:float):
        self.Thl=TempHotLeg
        self.CoolantVhl=28.3168
        self.Wc=18899.68208  
        self.CoolantDensity=732.18278
        self.Mhl=self.CoolantVhl*self.CoolantDensity
        self.time_consthl=self.Mhl/self.Wc

    def DThl(self,UpperPlenum:object):

        dtdThl=(UpperPlenum.Tup-self.Thl)/self.time_consthl
        return dtdThl

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


class ColdLeg():
    def __init__(self,TempColdLeg:float) -> None:
        self.Tcl=TempColdLeg
        self.Wc=18899.68208  
        self.CoolantDensity=732.18278
        self.CoolantVcl=56.63369
        self.Mcl=self.CoolantVcl*self.CoolantDensity
        self.time_constcl=self.Mcl/self.Wc

    def DTcl(self,TempRCPexit:float):
        """temperature exit will accessed through objects. but right now it's initilized as 
        a float data type """
        dtdTcl=(TempRCPexit-self.Tcl)/self.time_constcl
        return dtdTcl
    
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
    
class UpperPlenum():
    def __init__(self,TempUpperPlenum:float):
        self.Tup=TempUpperPlenum
        self.Wc=18899.68208  
        self.CoolantDensity=732.18278
        self.CoolantVup=38.9813
        self.Mup=self.CoolantDensity*self.CoolantVup
        self.time_constup=self.Mup/self.Wc

    def DTup(self,ReactorCoreInside:object):
        dtdTup=(ReactorCoreInside.Tmo6-self.Tup)/self.time_constup
        return dtdTup
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
        


class LowerPlenum():
    def __init__(self,TempLowerPlenum:float):
        self.Tlp=TempLowerPlenum
        self.Wc=18899.68208  
        self.CoolantDensity=732.18278
        self.CoolantVlp=50.7091
        self.Mlp=self.CoolantVlp*self.CoolantDensity
        self.time_constlp=self.Mlp/self.Wc

    def DTlp(self,ColdLeg:object):
        dtdTlp=(ColdLeg.Tcl-self.Tlp)/self.time_constlp
        return dtdTlp
    
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



ReactorCore=CoreInside(TemperatureFuel=[500,600,700],TemperatureModerator=[450,475,500,550,600,700],Power=1234,Precursor=.1)
t=0
dt=0.01
Time=[]
Temp=[]
while t<1000:

    ReactorCore.Tf1=ReactorCore.integrator(ReactorCore.DTf1,[],ReactorCore.Tf1,time_step=dt)
    ReactorCore.Tf2=ReactorCore.integrator(ReactorCore.DTf2,[],ReactorCore.Tf2,time_step=dt)
    ReactorCore.Tf3=ReactorCore.integrator(ReactorCore.DTf3,[],ReactorCore.Tf3,time_step=dt)
  
    print(ReactorCore.Tf1,' ',ReactorCore.Tf2,' ',ReactorCore.Tf3)

#Done 
