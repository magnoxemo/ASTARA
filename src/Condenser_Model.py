import numpy as np
from CoolProp.CoolProp import PropsSI


"""-------------------------------CONDNSER MODEL A-------------------------------"""


class Condenser_Model_A:
    def __init__(
        self,
        exhust_flow_rate_lp: float,
        exhust_temp_lp: float,
        exhust_temp_pressure: float,
        Condenser_outlet_temp: float,
        Temp_cold_water_in: float,
        Temp_hot_water_out: float,
    ):

        self.Wi = exhust_flow_rate_lp
        self.Ti = exhust_temp_lp
        self.time_const = 7
        self.Tout = Condenser_outlet_temp
        self.Rs = 0.4615e3
        self.Hot_well_diameter = 1
        self.volume = 3
        self.UA = 356.972e3
        self.Cp = 4.2e3

        self.hout = PropsSI("H", "T", self.Tout, "Q", 0, "water")
        self.hf = PropsSI("H", "T", self.Ti, "Q", 0, "water")
        self.hg = PropsSI("H", "T", self.Ti, "Q", 1, "water")
        self.hfg = self.hg - self.hf
        self.W2 = self.Wi * (self.hi - self.hf) / self.hfg
        self.W1 = self.Wi * (1 - ((self.hi - self.hf) / self.hfg))
        self.T_coldwater = Temp_cold_water_in
        self.T_hotwater = Temp_hot_water_out

        self.deltaT = (self.T_hotwater - self.T_coldwater) / np.log(
            (self.Ti - self.T_coldwater) / (self.Ti - self.T_hotwater)
        )
        self.hin = PropsSI("H", "T", self.Tin, "P", exhust_temp_pressure, "water")
        self.hcw = PropsSI("H", "T", self.T_coldwater, "Q", 0, "water")
        self.W3 = (self.UA * self.deltaT) / (self.hin - self.hcw)

    def DW_condensate(self):
        dtdWsc = (self.W2 - self.W3) / self.time_const
        return dtdWsc

    def DPs(self):

        dtdPs = self.Rs * (self.DWs() * self.Ti) / self.volume
        return dtdPs

    def Dhout(self):
        self.rho_w = PropsSI("D", "T", self.Tout, "Q", 0, "water")
        area = np.pi * self.Hot_well_diameter**2 / 4
        self.Lw = self.DW_condensate() / (self.rho_w * area)

        dtdho = (
            (self.W1 + self.W3) * (self.hf - self.hout) / (self.Lw * area * self.rho_w)
        )
        return dtdho

    def integrator(self, function, argsforfunction: [], intitial_cond, time_step):
        l = len(argsforfunction)

        if l == 0 or argsforfunction == None:
            return function() * time_step + intitial_cond
        elif l == 1:
            arg1 = argsforfunction[0]
            return function(arg1) * time_step + intitial_cond
        elif l == 2:
            arg1 = argsforfunction[0]
            arg2 = argsforfunction[1]
            return function(arg1, arg2) * time_step + intitial_cond
        elif l == 3:
            arg1 = argsforfunction[0]
            arg2 = argsforfunction[1]
            arg3 = argsforfunction[2]
            return function(arg1, arg2, arg3) * time_step + intitial_cond
        else:
            raise AttributeError(
                "agrs in your differential function were not correct! Fix them"
            )


"""-------------------------------CONDNSER MODEL B-------------------------------"""


class Condenser_Model_B:

    def __init__(
        self, steampressure: float, airpressure: float, hot_well_diametr: float
    ):

        self.volume = 3
        self.UA = 356.972e3
        self.Cp = 4.2e3
        self.Ps = steampressure
        self.Pa = airpressure
        self.Rs = 0.4615e3
        self.Ra = 0.287e3

        """ steam mass balance variables"""
        self.W_turbine = 4
        self.W_otherthanturbine = 0
        self.W_condensate = 4
        self.W_steamairout = 0

        self.T_steamin = 600
        self.T_steamout = 600

        """ cold water """
        self.W_coldwater = 107.881
        self.T_coldwater = 273 + 60
        self.Ww = 102310
        self.T_hotwater = 273 + 80
        self.hotwelldensity = 227
        self.Hs = PropsSI("H", "T", self.T_steamin, "P", self.Ps, "water")
        self.Hcw = PropsSI("H", "T", self.T_coldwater, "P", 101325, "water")
        self.deltaT = (self.T_hotwater - self.T_coldwater) / np.log(
            (self.T_steamin - self.T_coldwater) / (self.T_steamin - self.T_hotwater)
        )

        """ air zone """
        self.W_vaccumbreakvalve = 0
        self.W_air = 0
        self.W_steamgas = 0
        self.W_draincondenser = 0

        """ hot well water"""
        self.Hot_wellarea = np.pi * hot_well_diametr**2 / 4
        self.W_hotwell = 110
        self.W_bubblingoxygen = 100

        self.Wc = (self.UA * self.deltaT) / (self.Hs - self.Hcw)

    def DWs(self):

        self.R = (self.Pa * self.Rs) / (self.Pa * self.Ra + self.Ps * self.Rs)
        self.Wss = self.W_air * (1 - self.R)
        self.deltaT = (self.T_hotwater - self.T_coldwater) / np.log(
            (self.T_steamin - self.T_coldwater) / (self.T_steamin - self.T_hotwater)
        )
        self.Wc = (self.UA * self.deltaT) / (self.Hs - self.Hcw)
        dtdWsteam = self.W_turbine + self.W_otherthanturbine - self.Wc - self.Wss

        return dtdWsteam

    def DPs(self):
        dtdPs = self.Rs * (self.DWs() * self.T_steamin) / self.volume
        return dtdPs

    def Dhs(self, temp_otherthanturbine: float, pressure_otherthanturbine: float):

        self.Hs = PropsSI("H", "T", self.T_steamin, "P", self.Ps, "water")
        dtdWSHS = (
            self.W_turbine * self.Hs
            + self.W_otherthanturbine
            * PropsSI(
                "H", "T", temp_otherthanturbine, "P", pressure_otherthanturbine, "water"
            )
            - (self.Wc + self.Wss) * self.Hs
        ) / self.WaterLevel

        return dtdWSHS

    def DWa(self):
        dtdWa = (
            self.W_vaccumbreakvalve
            + self.W_draincondenser
            + self.W_steamgas
            - self.W_air
        )
        return dtdWa

    def DPa(self):
        dtdpa = self.DWa() * self.Ra * self.T_steamin / self.volume
        return dtdpa

    def TotalPressure(self):
        return self.Pa + self.Ps

    def Lw(self):
        return self.Ww / (self.Hot_wellarea * self.hotwelldensity)

    def DWw(self):
        dtdWhw = self.Wc + self.W_bubblingoxygen - self.Ww
        return dtdWhw

    def Dhw(self, Enthalpy_oxygen, pressure):
        dtdWwHw = (
            self.W_coldwater * self.Hcw
            + self.W_bubblingoxygen * Enthalpy_oxygen
            - self.Ww * PropsSI("H", self.T_steamout, "P", pressure, "water")
        )
        return dtdWwHw

    def DTtube(self):
        dtdTh = (
            self.UA * self.deltaT
            - self.W_coldwater * self.Cp * (self.T_hotwater - self.T_coldwater)
        ) / (self.W_coldwater * self.Cp)
        return dtdTh

    def integrator(self, function, argsforfunction: [], intitial_cond, time_step):
        l = len(argsforfunction)

        if l == 0 or argsforfunction == None:
            return function() * time_step + intitial_cond
        elif l == 1:
            arg1 = argsforfunction[0]
            return function(arg1) * time_step + intitial_cond
        elif l == 2:
            arg1 = argsforfunction[0]
            arg2 = argsforfunction[1]
            return function(arg1, arg2) * time_step + intitial_cond
        elif l == 3:
            arg1 = argsforfunction[0]
            arg2 = argsforfunction[1]
            arg3 = argsforfunction[2]
            return function(arg1, arg2, arg3) * time_step + intitial_cond
        else:
            raise AttributeError(
                "agrs in your differential function were not correct! Fix them"
            )


Condenser = Condenser_Model_B(
    steampressure=100e4, airpressure=89e3, hot_well_diametr=20
)

dt = 0.001
t = 0
i = 0
time = []
y = []

"""case studied: Sudden rain and temp drop  """
while t < 100:

    Condenser.Ps = Condenser.integrator(
        Condenser.DPs, [], intitial_cond=Condenser.Ps, time_step=dt
    )
    Condenser.W_turbine = Condenser.integrator(
        Condenser.DWs, [], intitial_cond=Condenser.W_turbine, time_step=dt
    )
    Condenser.T_hotwater = Condenser.integrator(
        Condenser.DTtube, [], intitial_cond=Condenser.T_hotwater, time_step=dt
    )

    t = t + dt
    i = i + 1

    if t > 50 and t < 800:
        Condenser.T_coldwater = Condenser.T_coldwater * np.exp(-(t - 50) / 100000)
    elif t > 800:
        Condenser.T_coldwater = (
            Condenser.T_coldwater * 0.001 * np.exp((t - 800) / 100000)
        )
    else:
        Condenser.T_coldwater = 290 + 60

    if i % 100 == 0:
        print("%3f" % t, " %.3e" % Condenser.Ps, " ", " %.3e" % Condenser.T_hotwater)
        time.append(t)
        y.append(Condenser.Ps)

import matplotlib.pyplot as plt

plt.plot(time, y)
plt.show()
