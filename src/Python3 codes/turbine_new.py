import numpy as np 
import scipy as sp 
from CoolProp.CoolProp import PropsSI

class Turbine_system():
    def __init__(self,steam_flow_rate_in:float,steam_generator_outlet_temp:float,rho_sd:float,
                 rho_hpex:float,T_hpex:float):
        self.Win=steam_flow_rate_in
        

        self.C1=0.2
        self.C2=1-self.C1
        self.Wsec=self.C2*self.Win
        self.Wmain=self.C1*self.Win

        """ Nozzle chest data """
        self.rho_c=38.8
        self.Pc=67.22e5
        self.hc=PropsSI("H",'D',self.rho_c,'P',self.Pc,'water')
        self.Tc=PropsSI('T','D',self.rho_c,'P',self.Pc,'water')
        self.Vc=212  #need to calculate 
        self.k1=0.02 #need to calculate

        '''  Reheater data   '''
        self.Tsd=steam_generator_outlet_temp
        self.rho_sd=rho_sd
        self.rho_r=2.4
        self.hr=2920.2e3
        self.Pr=PropsSI('P','H',self.hr,'D',self.rho_r,'water')

        ''' High pressure turbine'''
        self.rho_hpex=rho_hpex
        self.T_hpex=T_hpex
        self.hpex=PropsSI("H","D",self.rho_hpex,'T',self.T_hpex,'water')
        self.Cbhp=0.2            #this is at random. need no idea how to calculate it 
        self.Whp=self.Kch*np.sqrt(self.Pc*self.rho_c-self.Pr*self.rho_hpex)
        self.Whpex=299.2
        self.Tau1=2

        ''' Low pressure turbine '''

        self.Wlp=269.92
        self.Clp=0.02           #this is at random. need no idea how to calculate it 
        self.Wlpex=237.15
        self.hlpex=2457e3
        self.tau2=10

        ''' Feed Water heater High pressure '''
        self.hfwh1=425.4
        ''' feed water heater low pressure  '''



    def Throttle_valve(self):
        self.Wsec=self.C2*self.Win
        self.Wmain=self.C1*self.Win

    def Drho_c(self):
        self.Pc=PropsSI("P",'H',self.hc,'P',self.rho_c,'water')
        self.Wh1=self.Kch*np.sqrt(self.Pc*self.rho_c-self.Pr*self.rho_hpex)
        dtdrho_c=(self.Wmain-self.Wh1)/self.Vc
        return dtdrho_c 
    
    def W_mss(self):
        hf=PropsSI("H",'T',self.T_hpex,'Q',0,'water')
        hg=PropsSI("H",'T',self.T_hpex,'Q',0,'water')
        h_hpex=PropsSI("H",'T',self.T_hpex,'D',self.rho_hpex,'water')
        self.Wmss=self.Whpex*(h_hpex-hf)/(hg-hf)

    def W_msw(self):
        hf=PropsSI("H",'T',self.T_hpex,'Q',0,'water')
        hg=PropsSI("H",'T',self.T_hpex,'Q',0,'water')
        h_hpex=PropsSI("H",'T',self.T_hpex,'D',self.rho_hpex,'water')
        self.Wmsw=self.Whpex*(1-(h_hpex-hf)/(hg-hf))

    def Dh_c(self):
        hsd=PropsSI("H",'T',self.Tsd,'D',self.rho_sd,'water')
        self.Wh1=self.Kch*np.sqrt(self.Pc*self.rho_c-self.Pr*self.rho_hpex)
        a=(self.Wmain*hsd-self.whp1*self.hc)/(self.rho_c*self.Vc)
        b=(self.Pc*self.Drho_c())/self.rho_c**2
        dtdhc=(a+b)*(1/(1-self.k1))
        return dtdhc
    
    def DW_hpex(self):
        dtdWhpex=(self.Wh1*(1-self.Cbhp)-self.Whpex)/self.Tau1
        return dtdWhpex
    
    def DW_lpex(self):
        dtdWlpex=(self.Wlp*(1-self.Clp)-self.Wlpex)/self.Tau2
        return dtdWlpex
    

    def Torque(self):
        Torque_hp=self.Wh1*(self.hc- self.h_hpex)
        Torque_lp=self.Wlp*(self.hr-self.hlpex)

        return Torque_hp+Torque_lp
        
    def Dh_fwh1(self,hco:float,W_condenser_outlet:float):
        ''' Here needs some updates '''
        self.Wlp=self.Wr
        self.Wfw=W_condenser_outlet
        hh1=9232 # I don't know that value
        w1=self.Cbhp*self.Wh1+self.Wmss+self.Wr+self.Clp*self.Wlp
        H=(hco-self.hfwh1)
        Tau=60
        dtdh_fwh1=(w1*hh1/self.Wfw+H)/Tau
        return dtdh_fwh1
