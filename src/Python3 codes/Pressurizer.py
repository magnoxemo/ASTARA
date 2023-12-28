import numpy as np 
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d

class Pressurizer():
    def __init__(self,pressurizer_diameter:float,Pressure:float,Temp_coldleg:float):

        self.Pressure=Pressure
        self.diameter=pressurizer_diameter
        self.l=14.2524
        self.l_w=8.5527
        self.J_p=5.4027
        self.hw_bar=9.7209e5

        ''' m_spray and Q automatic variable and depends on the operation. I will left this on the user to denfine 
            when to operate them.
        '''
        self.m_spray=0
        self.Q=0
        self.area=self.diameter**2*np.pi/4
        '''
        -----------------------------------------------------------------------------------------
                                where the index j = 1 to N represent coolant nodes in the following
                                order, lower plenum, coolant node 1 and 2, upper plenum, hot-leg, inlet
                                plenum, PCL 1 and 2 and outlet plenum, and cold-leg. Based on the surge 
                                will be calculated
        -----------------------------------------------------------------------------------------------
        '''
        self.V_times_nu=np.array([.5991,0.1814,0.1814,1.3164,0.2752,0.277,0.022,0.6022,0.2776,0.1927])

        ''' Q--->  0  ----> represents water
            Q--->  1  ----> represents steam
        '''
        self.rho_w=PropsSI("D",'P',self.Pressure,'Q',0,'water')
        self.rho_s=PropsSI("D","P",self.Pressure,'Q',1,'water')
        self.Nu_w=1/PropsSI("D",'P',self.Pressure,'Q',0,'water')
        self.Nu_s=1/PropsSI("D",'P',self.Pressure,'Q',1,'water')
        self.h_w=PropsSI(("H",'P',self.Pressure,'Q',0,'water'))
        self.h_spr=PropsSI("H","T",Temp_coldleg,'Q',0,'water')

        '''------------------partial derivative co-efficient----------'''
        self.k1p=PropsSI('d(D)/d(T)|P','P',self.Pressure,'Q',0,'Water')
        self.k2p=PropsSI('d(D)/d(T)|P','P',self.Pressure,'Q',1,'Water')
        self.k3p=PropsSI('d(h)/d(T)|P','P',self.Pressure,'Q',1,'Water')
        
        '''-------------calculation of partital derivative ----------'''
        Pressure_=np.logspace(3,16.909,100000,base=np.e)
        density_=np.array(PropsSI("D",'P',Pressure_,'Q',0,'water'))
        nu_s=1/density_
        delta_nu_s=np.diff(nu_s)
        delta_pressure=np.diff(Pressure_)
        self.k4p=interp1d(delta_pressure,delta_nu_s)

        '''------------------------------------------------ '''

    def DPp(self,DTemp_list:list,Temp_coldleg):

        '''------------------surge calculation----------'''
        if len(DTemp_list)!=14:
            raise ValueError('You need all the temperature derivative in the core:')
        else:
            m_surge=np.sum(self.V_times_nu*np.array(DTemp_list))

        nu_s=1/PropsSI("D","P",self.Pressure,'Q',1,'water')
        h_spr=PropsSI("H","T",Temp_coldleg,'Q',0,'water')
        h_w=PropsSI(("H",'P',self.Pressure,'Q',0,'water'))
        rho_w=PropsSI("D",'P',self.Pressure,'Q',0,'water')
        rho_s=PropsSI("D","P",self.Pressure,'Q',1,'water')        

        C1p=(PropsSI("D",'P',self.Pressure,'Q',0,'water')/PropsSI("D","P",self.Pressure,'Q',1,'water'))-1
        
        k3p=PropsSI('d(h)/d(T)|P','P',self.Pressure,'Q',1,'Water')
        k4p=self.k4p(self.Pressure) #interpolated function 
        Vw=self.area*self.l_w
        Vs=self.area*(self.l-self.l_w)

        a=self.Q+m_surge*((self.Pressure*nu_s)/(self.J_p*self.C1p))
        b=self.m_spray*(h_spr-h_w+(self.hw_bar/C1p)+(self.Pressure/(rho_w*self.J_p*C1p)))
        c=Vw*rho_w*(k3p+(k4p*self.Pressure)/self.J_p)+(Vs*rho_s*k4p*self.Pressure-Vw)/self.J_p
        d=(self.hw_bar+self.Pressure*nu_s/self.J_p)

        dtdPp=(a+b)/(c+d)
        return dtdPp
    
    def Dlw(self,DTemp_list:list,Temp_coldleg):
        rho_s=PropsSI("D","P",self.Pressure,'Q',1,'water')

        '''------------------surge calculation----------'''
        if len(DTemp_list)!=14:
            raise ValueError('You need all the temperature derivative in the core:')
        else:
            m_surge=np.sum(self.V_times_nu*np.array(DTemp_list))

        k1p=PropsSI('d(D)/d(T)|P','P',self.Pressure,'Q',0,'Water')
        k2p=PropsSI('d(D)/d(T)|P','P',self.Pressure,'Q',1,'Water')
        C1p=(PropsSI("D",'P',self.Pressure,'Q',0,'water')/PropsSI("D","P",self.Pressure,'Q',1,'water'))-1
        C2p=self.area*(self.l-self.l_w)*(C1p+1)*k2p+self.area*self.l_w*k1p

        a=(self.area*(self.l-self.l_w)*k2p-C1p/C2p)*self.DPp(DTemp_list,Temp_coldleg)
        b=(m_surge/C1p)+((C2p*self.DPp(DTemp_list,Temp_coldleg))-m_surge-self.m_spray)/(C1p)**2
        dtdlw=(a+b)/rho_s*self.area

        return dtdlw
    
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


'''------------------------------------------------------------ Done & final ---------------------------------------------------- '''   
