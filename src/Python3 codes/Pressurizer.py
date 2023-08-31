import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d

class Pressurizer():

    def __init__(self,area:float,effective_length_of_pressurizer:float,water_level:float,pressure_at_pressurizer:float,
                 Temperature:float):
        
        """ 
        Wsp---> spray flow rate that will come from the any of the cold leg of any of the 4 steam generator 
        Wsr---> surge flow rate that will come from the any of the hot leg of any of the 4 steam generator 

        so when there is any annomalies we either control the Heater (Q), Wsp or The Wsr
        """
        
        #design parameters
        self.area=area
        self.L=effective_length_of_pressurizer
        self.Ppr=pressure_at_pressurizer
        self.Temp=Temperature

        #constants 

        pressure=np.linspace(700,2e7,num=20)
        pressureGrad=np.gradient(pressure)

        blank=[]

        for i in pressure:
            blank.append(PropsSI("V","P",i,'Q',0,'water'))

        self.Kp1=interp1d(pressure,np.gradient(blank)/pressureGrad)
        
        blank=[]

        for i in pressure:
            blank.append(PropsSI("V","P",i,'Q',1,'water'))

        self.Kp2=interp1d(pressure,np.gradient(blank)/pressureGrad)
        self.Kp3=interp1d(pressure,self.Kp2(pressure)-self.Kp1(pressure))
        
        blank=[]

        for i in pressure:
            blank.append(PropsSI("D","P",i,'Q',1,'water'))

        self.Kp4=interp1d(pressure,np.gradient(blank)/pressureGrad)
        
        blank=[]

        for i in pressure:
            blank.append(PropsSI("D","P",i,'Q',0,'water'))

        self.Kp5=interp1d(pressure,np.gradient(blank)/pressureGrad)

        blank=[]

        for i in pressure:
            blank.append(PropsSI("H","P",i,'Q',0,'water'))

        self.Kp6=interp1d(pressure,np.gradient(blank)/pressureGrad)
        
        blank=[]

        for i in pressure:
            blank.append(PropsSI("H","P",i,'Q',1,'water'))

        self.Kpg=interp1d(pressure,np.gradient(blank)/pressureGrad)
        self.Kp7=interp1d(pressure,self.Kpg(pressure)-self.Kp6(pressure))

        
        """ const[0]=d meu_f/dpr
            const[1]=d meu_g/dpr
            const[2]=d meu_fg/dpr
            const[3]=d rou_s/dpr
            const[4]=d rou_w/dpr
            const[5]=d h_f/dpr
            const[6]=d h_fg/dpr 
        """

        self.hf=PropsSI("H","P",self.Ppr,"Q",0,"water")
        self.hg=PropsSI("H","P",self.Ppr,"Q",1,"water")
        self.rouw=PropsSI("D","P",self.Ppr,"Q",0,"water")
        self.rous=PropsSI("D","P",self.Ppr,"Q",1,"water")

        self.hfg=self.hg-self.hf

        self.hsp=PropsSI("H","P",self.Ppr,'T',self.Temp,'water')

        self.vf=PropsSI("V","P",self.Ppr,"Q",0,"water")
        self.vg=PropsSI("V","P",self.Ppr,"Q",1,"water")
        self.vfg=self.vg-self.vf

        self.Wsp=0
        self.Wsr=0
        self.Q=0
        

        #variables

        self.Lw=water_level
        self.Vw=self.area*self.Lw

        #calculating some constant parameters 

        self.C1=self.rouw/self.rous -1
        self.C2=(self.area*(self.L-self.Lw))*(self.rouw/self.rous)*self.Kp5(self.Ppr) +self.Kp4(self.Ppr)*self.area*self.Lw
        

    
    def DLw(self):

        dtdLw=(1/(self.area*self.Lw))*(self.area*(self.L-self.Lw)*self.Kp5-self.C2/self.C1)*self.dPr()+self.Wsp/self.C1+self.Wco/self.C2
        self.Wco=(self.C2*self.DPr()-self.Wsr-self.Wsp)/self.C1
        return dtdLw

    def DPr(self):

        numerator=self.Q+self.Wsr*(self.Ppr*self.vg/(self.C1)+self.hfg/self.C1)\
        +self.Wsp*(self.hsp-self.hf+self.hfg/self.C1+self.Ppr*(self.vg+self.vfg)/(self.C1))

        denominator=self.Vw*self.rouw*((self.Kp6(self.Ppr))+self.Kp2(self.Ppr)*self.Ppr)+\
            (self.L-self.Lw)*self.area*self.rous*self.Kp2(self.Ppr)*self.Ppr-self.Vw +\
            (self.C2/self.C1)*(self.hfg+self.Ppr*(self.vf+self.vfg))

        dtdPr=numerator/denominator

        return dtdPr
    
    def PropertiesUpdater(self):

        self.hf=PropsSI("H","P",self.Ppr,"Q",0,"water")
        self.hg=PropsSI("H","P",self.Ppr,"Q",1,"water")
        self.rouw=PropsSI("D","P",self.Ppr,"Q",0,"water")
        self.rous=PropsSI("D","P",self.Ppr,"Q",1,"water")

        self.hfg=self.hg-self.hf

        self.hsp=PropsSI("H","P",self.Ppr,'T',self.Temp,'water')

        self.vf=PropsSI("V","P",self.Ppr,"Q",0,"water")
        self.vg=PropsSI("V","P",self.Ppr,"Q",1,"water")

        self.Wsp=0
        self.Wsr=0
        self.Vw=self.area*self.Lw

        #calculating some constant parameters 

        self.C1=self.rouw/self.rous -1
        self.C2=(self.area*(self.L-self.Lw))*(self.rouw/self.rous)*self.Kp5(self.Ppr) +self.Kp4(self.Ppr)*self.area*self.Lw

    def integrator(self,function,argsforfunction:None,intitial_cond,time_step):
        
        try:
            a=np.array(argsforfunction)
            l=len(a)
        except:
            pass

        if argsforfunction==None:
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
        
