import numpy as np

class PrimaryCoolantPump():

    def __init__(self,area_of_piping:float,piping_length:float,friction_factor:float,
                 moment_of_innertia:float,actual_pump_speed:float,norminal_flow_rate:float,actual_flow_rate:float,
                 power_delivered_to_shaft:float,const:list):

        self.area=area_of_piping
        self.L=piping_length
        self.K=friction_factor
        self.g=9.8
        self.water_density=700
        self.A=const[0]
        self.B=const[1]
        self.C=const[2]

        """ Initial condition """
        self.Np=actual_pump_speed
        self.Qp=actual_flow_rate
        self.Qs=norminal_flow_rate
        self.Pd=power_delivered_to_shaft

        self.beta=actual_flow_rate/norminal_flow_rate
        self.Hp=self.A*(self.Qp)**2+self.B*self.beta*self.Qp+self.C*self.beta**2
        self.I=moment_of_innertia


    def DevelopedPumpHead(self):

        self.Hp=self.A*(self.Qp)**2+self.B*self.beta*self.Qp+self.C*self.beta**2
    
    def DNp(self):

        dtdNp=(self.Pd-self.water_density*self.Qp*self.Hp)/(self.Np*self.I*4*np.pi**2)
        return dtdNp
    
    def DQ(self):
        dtdQ=(self.Qp-self.K*self.Qp**2)*self.g*self.area/self.L
        return dtdQ

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
