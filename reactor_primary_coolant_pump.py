import numpy as np


"""
Project NAME:
                         --------" ASTARA --a Nuclear Power Plant simulator "-------- 
PROGRAMMER:

EBNY WALID AHAMMED 
Undergrad Student (Level 4 term 1)
Dept of Nuclear Engineering 
University of Dhaka"""


class reactor_primary_coolant_pump:

    def __init__(
        self,
        area_of_piping: float,
        piping_length: float,
        friction_factor: float,
        # initial_contidions
        moment_of_innertia: float,
        norminal_pump_speed: float,
        actual_pump_speed: float,
        norminal_flow_rate: float,
        actual_flow_rate: float,
        power_delivered_to_shaft: float,
        torque_delivered_to_shaft: float,
        Total_hydrodynamic_head_loss: float,
        developed_pump_head: float,
        Hydrodynamic_torque: float,
        const: list,
    ):

        self.A = area_of_piping
        self.L = piping_length
        self.k = friction_factor
        self.g = 9.8
        self.water_density = 1000
        self.A = const[0]
        self.B = const[1]
        self.C = const[2]

        self.Np = actual_pump_speed
        self.Ns = norminal_pump_speed
        self.Qp = actual_flow_rate
        self.Qs = norminal_flow_rate
        self.Pd = power_delivered_to_shaft
        self.Td = torque_delivered_to_shaft
        self.Hd = developed_pump_head
        self.Ht = Total_hydrodynamic_head_loss
        self.I = moment_of_innertia
        self.Td = torque_delivered_to_shaft
        self.Th = Hydrodynamic_torque

        self.alpha = actual_flow_rate / norminal_flow_rate
        self.beta = actual_flow_rate / norminal_flow_rate

    def dqdt(self):

        self.Hp = (
            self.A * (self.Qp) ** 2
            + self.B * self.beta * self.Qp
            + self.C * self.beta**2
        )
        Dq = self.A * self.g * (self.Hp - self.Ht) / self.l

        return Dq

    def dNpdt(self):

        DNp = (self.Pd - self.water_density * self.Qp * self.Hp) / (
            self.Np * self.I * 4 * np.pi**2
        )
        return DNp

    def solver(self, function, y0, dt: float):

        f = function
        return y0 + function * dt
