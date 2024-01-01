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
        self.W_otherthanturbine=0
        self.W_condensate=4
        self.W_steamairout=0

        self.T_steamin=600
        self.T_steamout=400

        
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
        self.hotwelldensity=227
        self.Hs=PropsSI('H','T',self.T_steamin,'P',self.Ps,'water')
        self.Hcw=PropsSI('H','T',self.T_coldwater,'P',101325,'water')
        self.deltaT=(self.T_hotwater-self.T_coldwater)/np.log((self.T_steamin-self.T_coldwater)/(self.T_steamin-self.T_hotwater))
        print(self.deltaT)
        self.Wc=(self.UA*self.deltaT)/(self.Hs-self.Hcw)


    def DW_steam(self):

        self.R=(self.Pa*self.Rs)/(self.Pa*self.Ra+self.Ps*self.Rs)
        self.Wss=self.W_air*(1-self.R)
        dtdWsteam=self.W_turbine+self.W_otherthanturbine-self.W_condensate-self.Wss

        return dtdWsteam
    
    def DPs(self):
        dtdPs=self.Rs*(self.DW_steam()*self.T_steamin)/self.volume
        return dtdPs
    def Dsteamenthalpy(self,temp_otherthanturbine:float,pressure:float):
        dtdWSHS=self.W_turbine*self.Hs+self.W_otherthanturbine*PropsSI("H",'T',temp_otherthanturbine,'P',pressure,'water')-\
        (self.Wc+self.Wss)*self.Hs
        return dtdWSHS
    def DWa(self):
        dtdWa=self.W_vaccumbreakvalve+self.W_draincondenser+self.W_steamgas-self.W_air
        return dtdWa
    def TotalPressure(self):
        return self.Pa+self.Ps
    def DPa(self):
        dtdpa=self.DWa()*self.Ra*self.T_steamin/self.volume
        return dtdpa
    def WaterLevel(self):
        return self.W_hotwater/(self.Hot_wellarea*self.hotwelldensity)
    def DWhotwellwater(self):
        dtdWhw=self.W_condensate+self.W_bubblingoxygen-self.W_hotwater
    def DHotWellEnthalpy(self,Enthalpy_oxygen,pressure):
        dtdWwHw=self.W_coldwater*self.Hcw+self.W_bubblingoxygen*Enthalpy_oxygen-self.W_hotwater*PropsSI('H',self.T_steamout,'P',pressure,'water')
        return dtdWwHw
    def DT_hotwater(self):
        dtdTh=(self.UA*self.deltaT-self.W_coldwater*self.Cp*(self.T_hotwater-self.T_coldwater))/(self.W_coldwater*self.Cp)
        return dtdTh

    def integrator(self,function,argsforfunction:None,intitial_cond,time_step):
        l=len(argsforfunction)

        if l==0 or argsforfunction== None:
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

Condenser=Condenser(steampressure=10e4,airpressure=9e3)

dt=0.01
t=0
i=0
time=[]
y=[]

"""case studied: Sudden rain and temp drop  """
while t<1000:
    Condenser.Ps=Condenser.integrator(Condenser.DPs,[],intitial_cond=Condenser.Ps,time_step=dt)
    Condenser.W_turbine=Condenser.integrator(Condenser.DW_steam,[],intitial_cond=Condenser.W_turbine,time_step=dt)
    Condenser.T_hotwater=Condenser.integrator(Condenser.DT_hotwater,[],intitial_cond=Condenser.T_hotwater,time_step=dt)
    
    t=t+dt
    i=i+1
    if t>50 and t<800:
        Condenser.T_coldwater=Condenser.T_coldwater*np.exp(-(t-50)/1000000)
    elif t>800:
        Condenser.T_coldwater=Condenser.T_coldwater*0.01*np.exp((t-800)/100000)
    else:
        Condenser.T_coldwater=290+60
    if i%100==0:
        print ("%.3e" %Condenser.Ps,' ', "%.3e"%Condenser.T_hotwater)
        time.append(t)
        y.append(Condenser.T_hotwater)

import matplotlib.pyplot as plt
plt.plot(time,y)
plt.show()
            
