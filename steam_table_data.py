import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

class steam_data:
  
    def __init__(self):
        

        self.data_load()
        self.interpolate()


    def data_load(self):
    
        data=pd.read_csv('/home/ebny_walid/Downloads/H2O_TempSatda.csv')
        
        self.temp=np.array(data['  deg C'])+273
        self.p=np.array(data['     MPa'])*10**6
        self.vf=np.array(data['     vf'])*10**3
        self.vg=np.array(data['     vg'])*10**3
        self.uf=np.array(data['      uf'])*10**3
        self.ug=np.array(data['      ug'])*10**3
        self.hf=np.array(data['     hf'])*10**3
        self.hfg=np.array(data['     hfg      '])*10**3
        self.hfg=np.array(data['     hfg      '])*10**3

        self.sg=np.array(data['    sg'])*10**3
        self.sf=np.array(data['      sf'])*10**3
        self.sfg=np.array(data['   sfg'])*10**3

    def interpolate(self):
    

        self.entropy_fg=interp1d(self.temp,self.sfg)
        self.entropy_f=interp1d(self.temp,self.sf)
        self.enthalpy_fg=interp1d(self.temp,self.hfg)
        self.enthalpy_f=interp1d(self.temp,self.hf)
        self.internal_energy_g=interp1d(self.temp,self.ug)
        self.internal_energy_f=interp1d(self.temp,self.uf)
        self.vol_g=interp1d(self.temp,self.vg)
        self.vol_f=interp1d(self.temp,self.vf)
        self.pressure=interp1d(self.temp,self.p)

    def getPressure(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("Pressure data at this %.3f K Temperature range not found "  %Temperature)
        else:
            return self.pressure(Temperature)
    
    def getEntropy_f(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("Entropy_f data at this %.3f K Temperature range not found "  %Temperature)
        else:
            return self.entropy_f(Temperature)
    
    def getEntropy_fg(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("Entropy_fg at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.entropy_fg(Temperature)
    
    def getVolume_f(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("Volume_f data at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.vf(Temperature)
    
    def getVolume_g(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("Volume_g data at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.vg(Temperature)
    
    def getEnthalpy_f(self,Temperature:float):

        if Temperature<273 or Temperature>648:

            raise ValueError("Enthalpy_f data at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.enthalpy_f(Temperature)
    
    def getEnthalpy_fg(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("Enthalpy_fg data at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.enthalpy_fg(Temperature) 
    
    def getInternalEnergy_f(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("InternalEnergy_f data at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.internal_energy_f(Temperature)   
    
    def getInternalEnergy_g(self,Temperature:float):

        if Temperature<273 or Temperature>648:
            raise ValueError("InternalgEnergy_ data at this %.3f Temperature range not found "  %Temperature)
        else:
            return self.internal_energy_g(Temperature) 

Data_sheet=steam_data()
