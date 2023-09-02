import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
from CoolProp.CoolProp import PropsSI

def logo():
	print('          #                      #       ')
	print("         ####                  ####      ")
	print('        $$$$$$                $$$$$$     ')
	print('      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
	print('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('  %%%%%%%%@@@%%%%%%%%%%%%%%%%%%@@@%%%%%%%%%')
	print('%%%%%%%@@@@@@@@%%%%%%%%%%%%%%@@@@@@@@%%%%%%%%')
	print('%%%%%%%%@@@@@@%%%%%%%%%%%%%%%%@@@@@@%%%%%%%%%')
	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
	print('         %%%%%%%%%%%%%%%%%%%%%%%%%%%')
	print('             %%%%%%%%%%%%%%%%%%%%%')
	print('                 A   Z   O   G')
	print('    A Nuclear Power Plant Simulation code \n\n\n\n ')
	print('	   PROGRAMMER:')
	print("EBNY WALID AHAMMED ")
	print("Undergrad Student --Level 4 term 2")
	print("Dept of Nuclear Engineering")
	print("University of Dhaka")
		

class InletPlenum():
    
    def __init__(self,Temperature:float,pressure:float):
        """
        
        Mass     ---> stays in the pipe under steady state conditions
        Flow_rate---> hot_leg_flow_rate
        Theta    ---> temperature of the hot leg water 
        Pressure--> Reactor Pressure
        
        """
        self.Mass=10000
        self.Temperature=Temperature
        self.Pressure=pressure		 
        self.Flow_rate=4964.96
        self.density=PropsSI("D","T",self.Temperature,'P',self.Pressure,'water') 
        self.Theta=440          
        
    def DTpi(self):
        self.time_const=self.Flow_rate/self.Mass
        dtdTpi=(self.Theta-self.Temperature)*self.time_const
        
        """
        self.theta=reactor.T_hotleg
        self.theta should be coupled with the reactor 
        it will be reactor.T_hotleg
        
        """
        return dtdTpi
    
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



class PrimaryLump():

    def __init__(self,PrimaryLumpTemperature:list,MetalLumpTemperature:list,ProutTemperature:float,Pressure:float):
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

        self.Tpo=ProutTemperature
        self.Pressure=Pressure #Reactor Pressure or Pressure at the pressurizer

        self.Tavg=np.sum(PrimaryLumpTemperature)/len(PrimaryLumpTemperature)
        self.heat_capacity_1=PropsSI("C",'T',self.Tavg,'P',self.Pressure,'water')
        self.density=PropsSI("D",'T',self.Tavg,'P',self.Pressure,'water')

        ''' the heat capacity and the density needs to be constantly updated based on the temperature and the 
        pressure '''

        self.Primary_side_flim_conductance=25563
        self.Tube_metal_conductance=12263.68
        self.Tube_metal_conductance_subcool=11186.216
        self.Tube_metal_conductance_boiling=34068
        self.number_of_utube=3383
        self.inner_diameter=0.019685
        self.outer_diameter=0.022225
        self.length=10.831712
        self.frist_lump_length=1.05017116

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
    

class MetalLump():
    def __init__(self,PrimaryLumpTemperature:list,MetalLumpTemperature:list,Temperature_SFSL:float,Temperature_SFBL:float):
        
        self.heat_capacity_m=460.547802 #specific heat of the metal 
        self.number_of_utube=3383
        self.inner_diameter=0.019685
        self.outer_diameter=0.022225
        self.metaldensity= 8050.095
        self.length=10.831712
        self.frist_lump_length=1.05017116

        '''Mass calculation in the metal lump '''
        self.Mm=self.metaldensity*self.number_of_utube*self.length*np.pi*(self.outer_diameter**2-self.inner_diameter**2)/4
        self.Mm1=self.Mm*self.frist_lump_length/self.length
        self.Mm2=self.Mm*(self.length-self.frist_lump_length)/self.length
        self.Mm3=self.Mm2
        self.Mm4=self.Mm1

        '''      constant imported from the primary lump   '''

        self.Sm=np.pi*self.length*self.outer_diameter*self.number_of_utube
        self.Sm1=self.Sm*self.frist_lump_length/self.length
        self.Sm2=self.Sm*(self.length-self.frist_lump_length)/self.length
        self.Sm3=self.Sm2
        self.Sm4=self.Sm1

        self.Ums1=11186.216 #effective heat transfer co-efficient between water and steel 
        self.Ums2=14068
        self.Up1=25563

        

        '''initial conditions'''
        self.Td=Temperature_SFSL    #needs to be updated on a continiously 
        self.Tstat=Temperature_SFBL #needs to be updated on a continiously


        if len(PrimaryLumpTemperature)!=4:
            raise ValueError(" Initial condition error in metal lump")
        else:
            self.Tp1=PrimaryLumpTemperature[0]
            self.Tp2=PrimaryLumpTemperature[1]
            self.Tp3=PrimaryLumpTemperature[2]
            self.Tp4=PrimaryLumpTemperature[3]

        if len(MetalLumpTemperature)!=4:
            raise ValueError(" Initial condition error in metal lump")
        else:
            self.Tm1=MetalLumpTemperature[0]
            self.Tm2=MetalLumpTemperature[1]
            self.Tm3=MetalLumpTemperature[2]
            self.Tm4=MetalLumpTemperature[3]
    
    def DTm1(self,sub_cool_region:object):

        dtdTm1=self.Up1*self.Spm1*self.Tp1/(self.Mm1*self.heat_capacity_m)-\
        (self.Up1*self.Spm1+self.Ums1*self.Sm1)*self.Tm1/(self.Mm1*self.heat_capacity_m)\
        +self.Ums1*self.Sm1*(self.Td+self.Tstat)/(2*self.Mm1*self.heat_capacity_m)\
        -(self.Tm2-self.Tm1)*sub_cool_region.DLs1()/(2*self.frist_lump_length)

        return dtdTm1
    
    def DTm2(self,sub_cool_region:object):

        dtdTm2=self.Up1*self.Spm2*self.Tp2/(self.Mm2*self.heat_capacity_m)-(self.Up1*self.Spm2+self.Ums2*self.Sm2)*self.Tm2\
        /(self.Mm2*self.heat_capacity_m)+(self.Ums2*self.Sm2*self.Tstat)/(self.Mm2*self.heat_capacity_m)+(self.Tm2-self.Tm1)*\
        sub_cool_region.DLs1()/(2*(self.length-self.frist_lump_length))
        
        return dtdTm2
    
    def DTm3(self,sub_cool_region:object):

        dtdTm3=self.Up1*self.Spm2*self.Tp3/(self.Mm2*self.heat_capacity_m)-(self.Up1*self.Spm2+self.Ums2*self.Sm2)*self.Tm3\
        /(self.Mm2*self.heat_capacity_m)+(self.Ums2*self.Sm2*self.Tstat)/(self.Mm2*self.heat_capacity_m)+(self.Tm3-self.Tm4)*\
        sub_cool_region.DLs1()/(2*(self.length-self.frist_lump_length))
        
        return dtdTm3
    
    def DTm4(self,sub_cool_region:object):

        dtdTm4=self.Up1*self.Spm1*self.Tp4/(self.Mm1*self.heat_capacity_m)-\
        (self.Up1*self.Spm1+self.Ums1*self.Sm1)*self.Tm4/(self.Mm1*self.heat_capacity_m)\
        +self.Ums1*self.Sm1*(self.Td+self.Tstat)/(2*self.Mm1*self.heat_capacity_m)\
        -(self.Tm3-self.Tm4)*sub_cool_region.DLs1()/(2*self.frist_lump_length)

        return dtdTm4

    
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

class SubCooledRegion():
    def __init__(self,Tavg:float,MetalLump:object,PrimaryLump,HeaterConnectedToUTSG:object) :
        '''In and Out flow rate needs to be fixed '''

        self.area=5.63642501
        self.density=806.05092
        self.Cp2=4877.622

        "initial conditions "
        self.Ls1=1.05017116
        self.length=10.831712
        self.Tavg=Tavg

        self.W1=121
        self.W2=112
        self.MetalLump=MetalLump
        self.PrimaryLump=PrimaryLump
        self.HeaterConnectedToUTSG=HeaterConnectedToUTSG


    def DLs1(self):
        dtdLs1=(self.W1-self.W2)/(self.area*self.density)
        return dtdLs1
    
    def DTstat(self):

        k=self.MetalLump.Ums1*self.PrimaryLump.Pr2*self.Ls1*(self.MetalLump.Tm1+self.MetalLump.Tm4-self.MetalLump.Td-self.MetalLump.Tstat)+\
        self.W1*self.Cp2*self.MetalLump.Td-self.W2*self.Cp2*self.MetalLump.Tstat
        
        dtdTstat=(k/self.area*self.density)-(self.MetalLump.Td+self.MetalLump.Tstat)*self.DLs1()-self.Ls1*self.HeaterConnectedToUTSG.DTd()
        #DTd() will come from the heater 
        return dtdTstat
    
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
    
class BoilingRegion():
    def __init__(self,FlowRateOut:float,DowncomerTemp:float,BoilingTemp:float,SubCooledRegion:object):
        '''constants --> partial derivative const and enthalpies of
                                    hf.hg,hfg 
        '''
        self.SubCooledRegion=SubCooledRegion
        self.W3=FlowRateOut
        self.Td=DowncomerTemp
        self.Tstat=BoilingTemp
        self.Xe=0.2
        self.area=5.63642501
        
        self.hf=PropsSI('H','T',self.Tstat,'Q',0,'water')
        self.hg=PropsSI('H','T',self.Tstat,'Q',1,'water')
        self.hfg=self.hg-self.hf
        self.vf=PropsSI('V','T',self.Tstat,'Q',0,'water')
        self.vg=PropsSI('V','T',self.Tstat,'Q',0,'water')
        self.vfg=self.vg-self.vf
        self.density=PropsSI('D','T',self.Tstat,'Q',self.Xe,'water')

        """gradient constant Determination part"""
        Temp=np.linspace(300,600,num=200)
        Hf=[]
        Hg=[]

        for i in Temp:
            Hf.append(PropsSI("H","T",i,'Q',0,'water'))
            Hg.append(PropsSI("H","T",i,'Q',1,'water'))
                
        TempGrad=np.gradient(Temp)
        HfGrad=np.gradient(Hf)
        Hfg=np.array(Hg)-np.array(Hf)
        HfgGrad=np.gradient(Hfg)

        k1=HfGrad/TempGrad
        k2=HfgGrad/TempGrad
        self.dHfdTstatGrad=interp1d(Temp,k1)
        self.dHfgdTstatGrad=interp1d(Temp,k2)
        
    def DRoub(self):
        dtdRoub=((self.W1-self.W2)+self.density*self.area*self.SubCooledRegion.DLs1())/(self.SubCooledRegion.length-self.SubCooledRegion.Ls1)
        return dtdRoub
    
    def DXsteam(self,PrimaryLump:object,MetalLump:object,SubCooledRegion:object):

        k=MetalLump.Ums2*PrimaryLump.Pr2*(SubCooledRegion.length-SubCooledRegion.Ls1)*(MetalLump.Tm2+MetalLump.Tm3-2*SubCooledRegion.Tstat)+\
        SubCooledRegion.W2*self.hf-self.W3*(PropsSI('H','T',self.Tstat,'Q',0,'water')+self.Xe*(PropsSI('H','T',self.Tstat,'Q',1,'water')-\
        PropsSI('H','T',self.Tstat,'Q',0,'water')))

        self.hf=PropsSI('H','T',self.Tstat,'Q',0,'water')
        self.hg=PropsSI('H','T',self.Tstat,'Q',1,'water')

        dtdXe=((k/self.area)-self.density*(SubCooledRegion.length-SubCooledRegion.Ls1)*self.Xe*self.dHfgdTstatGrad(SubCooledRegion.Tstat)*SubCooledRegion.DTstat()/2+\
            self.density*(SubCooledRegion.length-SubCooledRegion.Ls1)*self.dHfdTstatGrad(SubCooledRegion.Tstat)*SubCooledRegion.DTstat()+\
            self.density*(SubCooledRegion.length-SubCooledRegion.Ls1)*(PropsSI('H','T',self.Tstat,'Q',self.Xe/2,'water'))*SubCooledRegion.DLs1()-\
            (SubCooledRegion.length-SubCooledRegion.Ls1)*(PropsSI('H','T',self.Tstat,'Q',self.Xe/2,'water'))*self.DRoub())/(self.density*(SubCooledRegion.length-SubCooledRegion.Ls1)*\
            (self.hg-self.hf))


        return dtdXe


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
    
class DrumRegion():
    def __init__(self,DrumWaterDTemperature:float,FeedWaterTemp:float,BoilingRegion:object,
                 PrimaryLump:object,MetalLump:object,SubCooledRegion:object):

        ''' Flow rates needs to be fixed '''
        self.area=10.2991

        #objects 
        self.BoilingRegion=BoilingRegion
        self.PrimaryLump=PrimaryLump
        self.MetalLump=MetalLump
        self.SubCooledRegion=SubCooledRegion

        """user defined value"""
        self.Xe=0.2        #steam quality
        self.water_level=2.935224 #Lw
        self.Wpi=4964.96   #hot_leg_flow_rate
        self.Wfi=470.226   #Trubine_outlet
        self.W1=2349.45    #SFSL
        self.W2=self.W1
        self.W3=self.W1
        self.W4=self.W1
        self.Wdw=1212
        #feed water is coming from the feed water pump after condensation 
        #so after this constructor it will be 
        #                                              "DrumRegion.Wfi=FeedWaterPump.outlet"

        Pressure=np.linspace(1e3,2.2e7,num=200)
        Vf=[]
        Vg=[]

        for i in Pressure:
            Vf.append(PropsSI("V","P",i,'Q',0,'water'))
            Vg.append(PropsSI("V","P",i,'Q',1,'water'))
                
        PressureGrad=np.gradient(Pressure)
        VfGrad=np.gradient(Vf)
        Vfg=np.array(Vf)-np.array(Vg)
        VfgGrad=np.gradient(Vfg)

        k1=VfGrad/PressureGrad
        k2=VfgGrad/PressureGrad
        self.dVfdPGrad=interp1d(Pressure,k1)
        self.dVfgdPGrad=interp1d(Pressure,k2)

        """initial conditions """
        self.Tw=DrumWaterDTemperature
        self.densityD=763.51
        self.Lw=2.935224 #Ldw
        self.Tfi=FeedWaterTemp

        """design parametrs of the DrumRegion """

        self.Vdr=124.28
        self.Pressure=5850053.972
        self.Cl=0.12232 #steam valve co efficient needs to be adjusted 

    def Dpressure(self):

        Vf=(PropsSI("V","P",self.Pressure,'Q',0,'water'))
        Vg=(PropsSI("V","P",self.Pressure,'Q',1,'water'))
        Vfg=Vf-Vg

        C1=-(self.dVfgdPGrad(self.Pressure)/(Vf+self.BoilingRegion.Xe*Vfg)**2)
        C2=-((self.dVfdPGrad(self.Pressure)+self.BoilingRegion.Xe*self.dVfgdPGrad(self.Pressure))/(Vf+self.BoilingRegion.Xe*Vfg)**2)

        dtdP=(((self.W2-self.W3)/self.Vdr)-C2*self.BoilingRegion.DXe(self.PrimaryLump,self.MetalLump,self.SubCooledRegion))/C1

        return dtdP
    
    def DLw(self):
        dtdlw=(-self.Wdw+(1-BoilingRegion.Xe)*self.W3-BoilingRegion.Xe*self.W3+self.Wfi)
        return dtdlw
    
    ''' done till here'''


    def DTw(self,MetalLump:object):

        val=self.Wfi*self.Tfi+(1-self.Xe)*self.W4*MetalLump.Tstat-self.W1*self.Tw 
        dtdTw=(val-self.density*self.area*self.Tw*self.DLw())/(self.density*self.area*self.water_level)

        return dtdTw
    
    def DDensityr(self,pressurechangerate,steamqualitychangerate):

        dtdrour=-(self.K1+self.K2*self.Xe)*(pressurechangerate)/(self.vf+self.Xe*self.vfg)**2\
        -self.vfg*steamqualitychangerate/((self.vf+self.Xe*self.vfg))**2
        return dtdrour 
    
    def DDensityg(self):

        dtdroug=(self.Xe*self.W4-self.Cl*self.Pressure+self.density*self.area*self.DLw())/(self.Vdr-self.area*self.water_level)
        return dtdroug

    
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
        
class DownComerRegion():
    def __init__(self,DownComerTemperature:float):  
        
        self.area=2.97376
        self.density=805.645
        self.Ld=10.8269
        self.Mass=self.area*self.Ld*self.density

        '''initial conditions '''


        self.Td=DownComerTemperature

    def DTd(self,DrumRegion:object):

        dtdTd=DrumRegion.W1*(DrumRegion.Tw-self.Td)/self.Mass
        return dtdTd
    
    
    def ConstOnUpdate(self,InletPlenum:object,PrimaryLump:object,MetalLump:object,SubcooledRegion:object,BoilingRegion:object,\
                            DrumRegion:object,DownCOmerRegion:object):
        
        MetalLump.Tstat=BoilingRegion.K5*DrumRegion.Pressure+771.205
        BoilingRegion.vf=BoilingRegion.K1*DrumRegion.Pressure+6453.186
        BoilingRegion.vfg=BoilingRegion.K2*DrumRegion.Pressure+0.288335
        BoilingRegion.density=1/(BoilingRegion.vf+BoilingRegion.vfg*BoilingRegion.Xe/2)
        DrumRegion.density=1/(BoilingRegion.vf+BoilingRegion.vfg*BoilingRegion.Xe)
        
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

    


InletPlenum=InletPlenum(430,pressure=10e6)
PrimaryLump=PrimaryLump(PrimaryLumpTemperature=[400,370,350,300],MetalLumpTemperature=[350,600,600,400],ProutTemperature=450,Pressure=10e6)

t=0
dt=0.001
T=[]
Temp=[]
Temp1=[]
dt=0.1

while t<100:

    InletPlenum.Temperature=InletPlenum.integrator(InletPlenum.DTpi,[],InletPlenum.Temperature,dt)
    PrimaryLump.Tp1=PrimaryLump.integrator(PrimaryLump.DTp1,argsforfunction=[InletPlenum],intitial_cond=PrimaryLump.Tp1,time_step=dt)

    Temp.append(PrimaryLump.Tp1-273)
    Temp1.append(InletPlenum.Temperature-273)
    T.append(t)
    t=dt+t
    #print("%.6f" %(InletPlenum.Temperature-273),"   ",'%.6f'%(PrimaryLump.Tp1-273))

from matplotlib import animation


def ani(i):
    plt.cla()
    plt.plot(T[:i],Temp[:i],color='red')
    plt.plot(T[:i],Temp1[:i],color='green')

ani = animation.FuncAnimation(plt.gcf(), ani,interval=1)

plt.show()
