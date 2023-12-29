import numpy as np 
import scipy as sp 
from CoolProp.CoolProp import PropsSI

class Reactor_Core():
    def __init__(self,TemperatureFuel:list,TemperatureModerator:list,Power:float,Precursor:float,Pressure:float,Temp_hotleg:float
                 ,Temp_coldleg:float,Temp_lowerplenum:float,Temp_upplerplenum:float):
        
        self.d=3.62712
        self.core_height=4.38912
        self.active_height=2
        self.CoolantDensity=PropsSI('D','P',Pressure,'Q',0,'water')
        self.rho_f=10920
        self.Cpc=PropsSI('C','P',Pressure,'Q',0,'water')   #coolant conductivity 
        self.area=14855.196091203/3
        self.U=1335                                         #needs more exect value conductivity

        #group const
        self.Beta1=0.000243
        self.Beta2=0.001363
        self.Beta3=0.001203
        self.Beta4=0.002605
        self.Beta5=0.000829
        self.Beta6=0.000166

        self.total_delayed_const=(self.Beta1+self.Beta2+self.Beta3+self.Beta4+self.Beta5+self.Beta6)
        
        #decay const
        self.Lamda1=0.0127
        self.Lamda2=0.0317
        self.Lamda3=0.115
        self.Lamda4=0.311
        self.Lamda5=1.40
        self.Lamda6=3.87


        self.NGT=20e-6
        self.Alpha_m=-1.8e-4
        self.Alpha_f=-2.16e-5
        self.Alpha_p=1.08e-4

        self.Fr=0.975                   #fission power factor 
        self.Cpf=.467e3                 #fuel heat conductivity 
        self.N=10693                    #number of fuel rods
        self.df=.95e-2                  #fuel diameter
        self.weight_of_a_fuel_rod=2
        
        #volumes 
        self.CoolantVcore=np.pi*(self.d/2)**2*self.core_height-np.pi*self.core_height*self.N*(self.df/2)**2

        self.FuelMass=self.N*self.rho_f*self.df**2*self.active_height/4
        self.Wc=1889.68208

        """initial conditions"""
        self.NominalPower=Power 
        self.Precursor=Precursor
        self.PowerRatio=1
        self.ExternalReactivity=0.00000 
        
         #control parameter will be used by control driving program 
         

        if len(TemperatureFuel)!=3:
            raise ValueError("Here should be three initial conditions!")
        else:
            self.Tf1=TemperatureFuel[0]
            self.Tf2=TemperatureFuel[1]
            self.Tf3=TemperatureFuel[2]
            self.Tf4=TemperatureFuel[3]
            self.Tf5=TemperatureFuel[4]
            
        if len(TemperatureModerator)!=6:
            raise ValueError('Here should be six initial conditions!')
        else:
            self.Tmo1=TemperatureModerator[0]
            self.Tmo2=TemperatureModerator[1]
            self.Tmo3=TemperatureModerator[2]
            self.Tmo4=TemperatureModerator[3]
            self.Tmo5=TemperatureModerator[4]
            self.Tmo6=TemperatureModerator[5]
            self.Tmo7=TemperatureModerator[6]
            self.Tmo8=TemperatureModerator[7]
            self.Tmo9=TemperatureModerator[8]
            self.Tmo10=TemperatureModerator[9]

        self.Tcl=Temp_coldleg
        self.Thl=Temp_hotleg
        self.Tup=Temp_upplerplenum
        self.Tlp=Temp_lowerplenum

        self.tempmod=np.sum(TemperatureModerator)
        self.tempfuel=np.sum(TemperatureFuel)
        self.Pressure=Pressure

        """ necessary calculations """
        """mass"""
        self.core_center_factor=0.5
        self.Mcl=1000
        self.Mlp=((self.core_height-self.active_height)*self.d**2*PropsSI("D","P",self.Pressure,"Q",0,'water'))*self.core_center_factor/4
        self.Mup=((self.core_height-self.active_height)*self.d**2*PropsSI("D","T",self.Tup,"Q",0,'water'))*(1-self.core_center_factor)/4
        self.Mhl=1000

        '''this values still not known. '''

        self.Mf=self.FuelMass/5                             #as there will  5  lumps 
        self.Mmo=self.CoolantVcore*self.CoolantDensity /10  #as there will 10  lumps

        #mass of the cold leg hot leg lower plenum

        """time const"""
        

        self.Lamda=self.total_delayed_const/(self.Beta1/self.Lamda1+self.Beta2/self.Lamda2+self.Beta3/self.Lamda3+self.Beta4/self.Lamda4\
                                             +self.Beta5/self.Lamda5+self.Beta6/self.Lamda6)
        
        self.reactivity=0

    def Reacivity(self,FuelTempSum:float,ModeratorTempSum:float,Pressure:float):
        
        self.reactivity=self.ExternalReactivity+self.Alpha_f*(FuelTempSum-self.tempfuel)/5+self.Alpha_m*(ModeratorTempSum-self.tempmod)/10+self.Alpha_p(Pressure-self.Pressure)
    
    def DPower(self):
        dtdPP0=((self.reactivity-self.total_delayed_const)*self.PowerRatio)/self.NGT+self.Lamda*self.Precursor
        return dtdPP0
        
    def DPrecoursor(self):
        dtdc=self.total_delayed_const*self.PowerRatio/self.NGT-self.Lamda*self.Precursor
        return dtdc
        
    def DTf1(self):
        dtdTf1=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo1-self.Tf1)/(self.Mf*self.Cpf)
        return dtdTf1
        
    def DTf2(self):
        dtdTf2=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.U*self.area*(self.Tmo3-self.Tf2)/(self.Mf*self.Cpf)
        return dtdTf2
        
    def DTf3(self):
        dtdTf3=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo5-self.Tf3)/(self.Mf*self.Cpf)
        return dtdTf3
    
    def DTf4(self):
        dtdTf2=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo7-self.Tf4)/(self.Mf*self.Cpf)
        return dtdTf2
        
    def DTf5(self):
        dtdTf3=self.Fr*self.NominalPower*self.PowerRatio/(self.Mf*self.Cpf)+self.h*self.area*(self.Tmo9-self.Tf5)/(self.Mf*self.Cpf)
        return dtdTf3
    
    def DTmo1(self,LowerPlenum:object):

        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.U*self.area*(self.Tf1-self.Tmo1)/(self.Mmo*self.Cpc)+(LowerPlenum.Tlp-self.Tmo1)/self.time_constmo
        dtdTm01=a+b
        return dtdTm01 
        
    def DTmo2(self):

        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.U*self.area*(self.Tf1-self.Tmo1)/(self.Mmo*self.Cpc)+(self.Tmo1-self.Tmo2)/self.time_constmo

        dtdTmo2=a+b

        return dtdTmo2 
    

    def DTmo3(self):
        
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.U*self.area*(self.Tf2-self.Tmo3)/(self.Mmo*self.Cpc)+(self.Tmo2-self.Tmo3)/self.time_constmo
        dtdTm03=a+b      
        return dtdTm03
    
        
    def DTmo4(self):
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.h*self.area*(self.Tf2-self.Tmo3)/(self.Mmo*self.Cpc)+(self.Tmo3-self.Tmo4)/self.time_constmo
        dtdTm04=a+b
        return dtdTm04
    
    def DTmo5(self):
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.h*self.area*(self.Tf3-self.Tmo5)/(self.Mmo*self.Cpc)+(self.Tmo4-self.Tmo5)/self.time_constmo
        dtdTm05=a+b
        return dtdTm05
    
    def DTmo6(self):

        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.h*self.area*(self.Tf3-self.Tmo5)/(self.Mmo*self.Cpc)+(self.Tmo5-self.Tmo6)/self.time_constmo

        dtdTmo2=a+b

        return dtdTmo2 
    

    def DTmo7(self):
        
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=(self.Tmo6-self.Tmo7)/self.time_constmo+self.h*self.area*(self.Tf4-self.Tmo7)/(self.Mmo*self.Cpc)
        dtdTm03=a+b      
        return dtdTm03
    
        
    def DTmo8(self):
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.U*self.area*(self.Tf4-self.Tmo7)/(self.Mmo*self.Cpc)+(self.Tmo7-self.Tmo8)/self.time_constmo
        dtdTm04=a+b
        return dtdTm04
    
    def DTmo9(self):
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.U*self.area*(self.Tf5-self.Tmo9)/(self.Mmo*self.Cpc)+(self.Tmo8-self.Tmo9)/self.time_constmo
        dtdTm05=a+b
        return dtdTm05
        
    def DTmo10(self):
        self.time_constmo=self.Mmo/(self.Wc*2)
        a=(1-self.Fr)*self.NominalPower*self.PowerRatio/(self.Mmo*self.Cpc)
        b=self.h*self.area*(self.Tf5-self.Tmo9)/(self.Mmo*self.Cpc)+(self.Tmo9-self.Tmo10)/self.time_constmo
        dtdTm06=a+b
        return dtdTm06
    
    def DTcl(self,Temp_RCP:float):

        dtdTcl=self.Wc*(Temp_RCP-self.Tcl)/self.Mcl
        return dtdTcl
   
    def DThl(self,):

        dtdTcl=self.Wc*(self.Tup-self.Thl)/self.Mhl
        return dtdTcl
    
    def DTup(self):
        self.Mup=((self.core_height-self.active_height)*self.d**2*PropsSI("D","T",self.Tup,"Q",0,'water'))*(1-self.core_center_factor)/4
        dtdTcl=self.Wc*(self.Tmo10-self.Tup)/self.Mup
        return dtdTcl
    
    def DTlp(self):
        self.Mlp=((self.core_height-self.active_height)*self.d**2*PropsSI("D","P",self.Pressure,"Q",0,'water'))*self.core_center_factor/4
        dtdTcl=self.Wc*(self.Tcl-self.Tlp)/self.Mlp
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


'''------------------------------------------------------------ Done---------------------------------------------------- ''' 



    
