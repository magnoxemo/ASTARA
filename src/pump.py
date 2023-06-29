import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sp 
from Reactor_MODEL import Reactor
from UTSG_MODEL import UTSG


class reactor_primary_coolant_pump():

    def __init__(self,area_of_piping:float,piping_length:float,friction_factor:float,
                 #initial_contidions
                 moment_of_innertia:float,norminal_pump_speed:float,actual_pump_speed:float,norminal_flow_rate:float,actual_flow_rate:float,
                 power_delivered_to_shaft:float,torque_delivered_to_shaft:float,Total_hydrodynamic_head_loss:float,
                 developed_pump_head:float,Torque_delivered_to_pump_shaft:float,Hydrodynamic_torque:float):

        self.A=area_of_piping
        self.L=piping_length
        self.k=friction_factor
        self.g=9.8
        self.water_density=1000


        self.Np=actual_pump_speed
        self.Ns=norminal_pump_speed
        self.Qp=actual_flow_rate
        self.Qs=norminal_flow_rate
        self.Pd=power_delivered_to_shaft
        self.Td=torque_delivered_to_shaft
        self.Hd=developed_pump_head
        self.Hl=Total_hydrodynamic_head_loss
        self.I=moment_of_innertia
        self.Td=torque_delivered_to_shaft
        self.Th=Hydrodynamic_torque

