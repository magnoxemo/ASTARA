import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sp 
from Reactor_MODEL import Reactor
from UTSG_MODEL import UTSG

class pump():

    def __init__(self,nominal_speed:float,
                 actual_speed:float,friction_factor:float,
                 primar_loop_pipe_length:float,
                 effective_area_piping_side:float) -> None:

        self.density_water=1000
        self.friction_const=friction_factor
        self.gravitaion=9.8
        self.nominal_pump_speed=nominal_speed
        self.actual_pump_speed=actual_speed
        self.L=primar_loop_pipe_length
        self.A_eff=effective_area_piping_side

        pass
