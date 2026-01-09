import numpy as np
import matplotlib.pyplot as plt
from CoolProp import CoolProp
form Base import Base

class reactor_primary_coolant_pump(Base):

    def __init__(
        self,
        area_of_piping: float,
        piping_length: float,
        friction_factor: float,
        moment_of_innertia: float,
        norminal_pump_speed: float,
        actual_pump_speed: float,
        norminal_flow_rate: float,
        actual_flow_rate: float,
        power_delivered_to_shaft: float,
        developed_pump_head: float,
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
        self.Hd = developed_pump_head
        self.I = moment_of_innertia

    def DQdt(self):

        beta = self.Np / self.Ns
        self.Hp = self.A * (self.Qp) ** 2 + self.B * beta * self.Qp + self.C * beta**2
        Dq = self.A * self.g * (self.Hp - self.k * self.Qp**2) / self.L

        return Dq

    def DNpdt(self):
        self.Td = self.Pd / (np.pi * 2 * self.Np)
        self.Th = (
            self.water_density * self.g * self.Qp * self.Hd / (np.pi * 2 * self.Np)
        )
        DNp = self.Td - self.Th
        return DNp


reactor_primary_coolant_pump = reactor_primary_coolant_pump(
    area_of_piping=10,
    piping_length=13,
    friction_factor=1,
    moment_of_innertia=800,
    actual_flow_rate=850,
    norminal_flow_rate=8498,
    norminal_pump_speed=1015450,
    actual_pump_speed=842,
    power_delivered_to_shaft=6,
    developed_pump_head=282,
    const=[0.001, -0.02, 0.001],
)


dt = 0.01
t = 0

T = []
n = []
q = []
while t < 100:

    Q = reactor_primary_coolant_pump.integrator(
        reactor_primary_coolant_pump.DQdt,
        [],
        intitial_cond=reactor_primary_coolant_pump.Qp,
        time_step=dt,
    )
    N = reactor_primary_coolant_pump.integrator(
        reactor_primary_coolant_pump.DNpdt,
        [],
        intitial_cond=reactor_primary_coolant_pump.Np,
        time_step=dt,
    )
    reactor_primary_coolant_pump.Np = N
    reactor_primary_coolant_pump.Qp = Q
    T.append(t)
    n.append(N)
    q.append(Q)
    if t > 20 and t < 60:
        reactor_primary_coolant_pump.L = 12
        reactor_primary_coolant_pump.Np = 3000
    else:
        reactor_primary_coolant_pump.Np = 3099
        reactor_primary_coolant_pump.L = 900
    print(Q, "", N)
    t = t + dt


plt.plot(T, q)
plt.plot(T, n)
plt.show()
