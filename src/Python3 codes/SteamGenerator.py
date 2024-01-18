import numpy as np 
import scipy as sp 
from CoolProp.CoolProp import PropsSI

"""    -------------------------- Model Begin -----------------------------------------"""

class u_tube_steam_generator(): 
    def __init__(self,primary_coolant_inlet_temperature:float,primary_coolant_outlet_temperature:float,feed_water_inlet_temperature:float,drum_water_temp:float,
                 avg_sub_cool_temp:float,feed_water_flow_rate:float,PrimaryLumpTemperature:list,MetalLumpTemperature:list,
                 Reactor_Pressure:float,Steam_pressure:float):
        
        '''------------------------       design parameters        -------------------'''
        self.N=3388
        self.L=10.83
        self.L_w=3.057   #sub cool region height
        self.Ldw=10.83   #Drum water level height
        self.R_in=0.019685 #needs to be checked 
        self.R_out=0.022225  #needs to be checked
        self.k=15
        self.Pressure_r=Reactor_Pressure
        self.Pressure_s=Steam_pressure
        self.Cl=80e-6
        '''------------------------        area and volume          ------------------'''
        self.P_r1=2*np.pi*self.R_in
        self.P_r2=2*np.pi*self.R_out
        self.Ap=self.N*np.pi*self.R_in**2
        self.Afs=5.63643
        self.Ad=9.39528444
        self.rho_m=8050                #needs the confirmation 
        self.rho_b=PropsSI('D',"P",self.Pressure_s,'Q',0.99/2,'water')
        self.Vp=self.Ap*self.L
        self.Vr=13.2523
        self.Vdr=124.55748301
        self.Vpi=(0.5*self.Vp-self.Ap*self.L)

        '''------------------- coolant and metal conductivity ------------------------'''
        if len(PrimaryLumpTemperature)!=4:
            raise ValueError("Initial condition error!")
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

        self.Tfi=feed_water_inlet_temperature
        self.Tpi=primary_coolant_inlet_temperature
        self.Tpo=primary_coolant_outlet_temperature
        self.Ts1=avg_sub_cool_temp
        self.Tsat=545.31667
        self.Td=535.76111
        self.Tdw=drum_water_temp

        self.Cp1=PropsSI("C","T",self.Tp1,'Q',0,'water') #reactor pressure instead of the Q= 0
        self.Cm=460

        '''-------------------------       flow rates             ----------------------'''
        self.Win=4950
        self.Wp1=self.Win
        self.Wpout=self.Wp1
        '''----------------------------------secondary loop------------------------------'''
        self.W1=892         #not sure! 
        self.W2=self.W1
        self.W3=self.W1
        self.W4=self.W1
        self.Wfi=feed_water_flow_rate
        self.rho_r=7.94*16.01844634
        '''-------------------------      mass calculation        ----------------------'''
        self.mpi=PropsSI('D','T',primary_coolant_inlet_temperature,'P',self.Pressure_r,'water')*self.Vpi #reactor pressure instead of the Q= 0
        self.mp1=PropsSI('D','T',self.Tp1,'P',self.Pressure_r,'water')*self.P_r1*self.L_w #reactor pressure instead of the Q= 0
        self.mp2=PropsSI('D','T',self.Tp2,'P',self.Pressure_r,'water')*self.P_r1*(self.L-self.L_w) #reactor pressure instead of the Q= 0
        self.mp3=PropsSI('D','T',self.Tp3,'P',self.Pressure_r,'water')*self.P_r1*(self.L-self.L_w)
        self.mp4=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')*self.P_r1*(self.L_w)

        self.mm1=self.N*np.pi*self.L_w*(self.R_out**2-self.R_in**2)*self.rho_m
        self.mm2=self.N*np.pi*(self.L-self.L_w)*(self.R_out**2-self.R_in**2)*self.rho_m

        self.hi=7.098
        self.hd=4.974
        self.hb=10.618

        self.Upm=1/((1/self.hi)+(self.R_in/self.k)*(np.log10((self.R_out+self.R_in)/((self.R_in*2)))))
        self.Ums1=1/((1/self.hd)+(self.R_out/self.k)*(np.log10(2*self.R_out/(self.R_out+self.R_in))))
        self.Ums2=1/((1/self.hb)+(self.R_out/self.k)*(np.log10(2*self.R_out/(self.R_out+self.R_in))))

    
    '''--------------------------------------------primary loop begin -------------------------------'''
  
    def DTpi(self,Hot_leg_temp:float):
        
        dtdTpi=self.Win*(Hot_leg_temp-self.Tpi)/self.mpi
        return dtdTpi
    
    def DTpo(self):
        dtdTpo=(self.Tp4-self.Tpo)*self.Win/(self.mpi)
        return dtdTpo
    
    def DLs1(self):
        rho_p=PropsSI('D','T',self.Tp1,'Q',0,'water')
        dtdLs1=(self.Win-self.Wp1)/(rho_p*self.Ap)
        return dtdLs1

    def DTp1(self):

        rho_p=PropsSI('D','T',self.Tp1,'P',self.Pressure_r,'water')
        Cp1=PropsSI('C','T',self.Tp1,'P',self.Pressure_r,'water')
        mp1=PropsSI('D','T',self.Tp1,'P',self.Pressure_r,'water')*self.P_r1*(self.L_w)

        a=self.Win*(self.Tpi-self.Tp1)/(rho_p*self.Ap*self.L_w)
        b=self.Upm*self.P_r1*self.L_w*self.N*(self.Tm1-self.Tp1)
        c=mp1*Cp1
        
        dtdTp1=a+b/c
        return dtdTp1
    
    def DTp2(self):

        rho_p=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')
        Cp1=PropsSI('C','T',self.Tp4,'P',self.Pressure_r,'water')
        mp2=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')*self.P_r1*(self.L-self.L_w)



        a=self.Win*(self.Tp1-self.Tp2)/mp2
        b=-(self.Upm*self.P_r1*self.N*(self.L-self.L_w)*self.N*(self.Tp2-self.Tm2))/(mp2*Cp1)
        c=-((self.Tp1-self.Tp2))/(self.L-self.L_w)
        d=self.DLs1()

        dtdTp2=a+b+c*d

        return dtdTp2
    
    def DTp3(self):

        rho_p=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')
        Cp1=PropsSI('C','T',self.Tp4,'P',self.Pressure_r,'water')
        mp3=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')*self.P_r1*(self.L-self.L_w)

        a=self.Win*(self.Tp2-self.Tp3)/(rho_p*self.Ap*(self.L-self.L_w))
        b=self.Upm*self.P_r1*self.N*(self.L-self.L_w)*(self.Tm3-self.Tp3)
        c=mp3*Cp1

        dtdTp3=a+b/c
        return dtdTp3   
          
    
    def DTp4(self):

        rho_p=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')
        Cp1=PropsSI('C','T',self.Tp4,'P',self.Pressure_r,'water')
        mp4=PropsSI('D','T',self.Tp4,'P',self.Pressure_r,'water')*self.P_r1*(self.L_w)

        a=self.Win*(self.Tp3-self.Tp4)/(rho_p*self.Ap*self.L_w)
        b=self.Upm*self.P_r1*self.N*(self.L_w)*(self.Tm4-self.Tp4)/(mp4*Cp1)
        c=(self.Tp3-self.Tp4)/self.L_w
        dtdTp4=a+b+c*self.DLs1()
        return dtdTp4
    
    '''--------------------------------------------primary loop done -------------------------------'''
    '''---------------------------------------------------------------------------------------------'''
    

    '''-------------------------------------------- Metal lump begin -------------------------------'''
    
    def DTm1(self):

        a=(self.Upm*self.N*self.P_r1*self.L_w*(self.Tp1-self.Tm1)-self.Ums1*self.N*self.P_r2*self.L_w*(self.Tm1-(self.Td+self.Tsat)/2))/(self.mm1*self.Cm)
        b=-(self.Tm1-self.Tm2)*0.5/self.L_w
        c=self.DLs1()
        dtdTm1=a+b*c
        return dtdTm1
    
    def DTm2(self):

        a=self.Upm*self.N*self.P_r1*(self.L-self.L_w)*self.Tp2/(self.mm2*self.Cm)
        b=-(self.Upm*self.N*self.P_r1*(self.L-self.L_w)+self.Ums2*self.N*self.P_r2*(self.L-self.L_w))*self.Tm2/(self.mm2*self.Cm)
        c=self.Ums2*self.N*self.P_r2*(self.L-self.L_w)*self.Tsat/(self.mm2*self.Cm)
        d=(self.Tm2-self.Tm1)/(2*(self.L-self.L_w))

        dtdTm2=a+b+c+d*self.DLs1()
        
        return dtdTm2
    
    def DTm3(self):

        a=self.Upm*self.N*self.P_r1*(self.L-self.L_w)*self.Tp2/(self.mm2*self.Cm)
        b=-(self.Upm*self.N*self.P_r1*(self.L-self.L_w)+self.Ums2*self.N*self.P_r2*(self.L-self.L_w))*self.Tm3/(self.mm2*self.Cm)
        c=self.Ums2*self.N*self.P_r2*(self.L-self.L_w)*self.Tsat/(self.mm2*self.Cm)
        d=(self.Tm3-self.Tm4)/(2*(self.L-self.L_w))

        dtdTm3=a+b+c+d*self.DLs1()
        
        return dtdTm3
    
    def DTm4(self):

        a=self.Upm*self.N*self.P_r1*(self.L_w)*self.Tp4/(self.mm1*self.Cm)
        b=-(self.Upm*self.N*self.P_r1*(self.L_w)+self.Ums1*self.N*self.P_r2*(self.L_w))*self.Tm4/(self.mm1*self.Cm)
        c=-self.Ums1*self.N*self.P_r2*(self.L_w)*self.Ts1/(self.mm1*self.Cm)
        d=((self.Tm3-self.Tm4)/(2*(self.L_w)))*self.DLs1()

        dtdTm4=a+b+c+d
        return dtdTm4
    

    '''-------------------------------------------- Metal lump done -------------------------------'''

    ''' ------------------------------------------------------------------------------------------'''
    ''' ----------------------------------------- secondary lump-----------------------------------'''
    
    def DTs1(self):

        Cp2=PropsSI("C","T",self.Ts1,'Q',0,'water')
        rho_s1=PropsSI("D",'T',self.Td,'Q',0,'water')
        self.W2=self.W1-self.DLs1()*rho_s1*self.Ap
        a=self.Ums1*self.N*self.P_r2*self.L_w*(self.Tm1+self.Tm4-2*self.Ts1)+self.W1*Cp2*self.Td-self.W2*Cp2*self.Tsat
        b=self.Afs*rho_s1*self.L_w*Cp2

        dtdTs1=a/b
        return dtdTs1

    def Drho_b(self):

        a=(self.W2-self.W3)/(self.Afs*(self.L-self.L_w))
        b=self.rho_b*self.DLs1()/(self.L-self.L_w)

        dtdrho_b=a+b
        return dtdrho_b
    
    def Dh_b(self):
        a=self.Ums2*self.N*self.P_r2*(self.L-self.L_w)*(self.Tm2-self.Tsat)+self.Ums2*self.N*self.P_r2*(self.L-self.L_w)*(self.Tm3-self.Tsat)
        b=self.W2*PropsSI('H','T',self.Ts1,'Q',0,'water')-self.W3*PropsSI('H','T',self.Ts1,'Q',0.233,'water')
        c=self.rho_b*self.Afs*self.hb*self.DLs1()
        f=(self.L-self.L_w)*self.hb*self.Drho_b()*self.Afs
        d=self.rho_b*self.Afs*(self.L-self.L_w)

        dtdh_b=(a+b+c-f)/d
        return dtdh_b
    
    def Drho_r(self):
        dtdrho_r=(self.W3-self.W2)/self.Vr
        return dtdrho_r
    
    def DLdw(self):
        Xe=0.2
        a=-self.W1-(1-Xe)*self.W3+self.Wfi
        rho_d=PropsSI('D',"T",self.Tfi,'Q',0,'water')
        dtdLdw=a/(rho_d*self.Ad)
        return dtdLdw
    
    def DTdw(self):
        Xe=0.2
        self.rho_g=PropsSI("D",' T',self.Tdw,"Q",Xe,'water')
        rho_d=PropsSI('D',"T",self.Tfi,'Q',0,'water')
        dtdTdw=(self.Wfi*self.Tfi-(1-Xe)*self.W4*self.Tsat-self.W1*self.Tdw-rho_d*self.Ad*self.DLdw())/(self.Ad*rho_d*self.Ldw)
        return dtdTdw
    
    def Drho_g(self):
        Xe=0.2
        self.rho_g=PropsSI("D",' T',self.Tdw,"Q",Xe,'water')
        a=self.W4*Xe-self.Cl*self.Pressure_s+self.rho_g*self.Ad*self.DLdw()
        b=self.Vdr-self.Ad*self.Ldw 
        dtdrho_g=a/b
        return dtdrho_g
    
    '''------------------------------- Down Comer region-------------------------'''

    def DTd(self):
        rho_d=PropsSI('D',"T",self.Tfi,'Q',0,'water')
        dtdTd=(self.Tdw-self.Td)*self.W1/(self.Ldw*self.Ad*rho_d)
        return dtdTd

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
