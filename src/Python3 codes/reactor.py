import numpy as np
import scipy as sp
from CoolProp.CoolProp import PropsSI


class Reactor_Core:
    def __init__(
        self,
        TemperatureFuel: list,
        TemperatureModerator: list,
        Power: float,
        Precursor: float,
        Pressure: float,
        Temp_hotleg: float,
        Temp_coldleg: float,
        Temp_lowerplenum: float,
        Temp_upplerplenum: float,
    ):
        self.d = 3.62712
        self.core_height = 4.38912
        self.CoolantDensity = PropsSI("D", "P", Pressure, "Q", 0, "water")
        self.active_height = 2
        self.df = 0.95e-2
        self.gas_gap = 0.057e-2
        self.N_rods = 10693
        self.N_fue_rods = 9768
        self.rho_f = 10960
        # self.area=np.pi*self.N*self.active_height*self.df**2/4
        self.U = 1370  # needs more exect value conductivity

        # group const
        self.Beta1 = 0.000215
        self.Beta2 = 0.001424
        self.Beta3 = 0.001274
        self.Beta4 = 0.002568
        self.Beta5 = 7.483e-4
        self.Beta6 = 2.73e-4

        self.total_delayed_const = (
            self.Beta1 + self.Beta2 + self.Beta3 + self.Beta4 + self.Beta5 + self.Beta6
        )

        # decay const

        self.Lamda1 = 1.2437e-2
        self.Lamda2 = 3.05e-2
        self.Lamda3 = 0.1141
        self.Lamda4 = 0.3013
        self.Lamda5 = 1.12866
        self.Lamda6 = 3.0130

        self.NGT = 2 * 10**-5
        self.Alpha_m = -1.1192 * 10**-7
        self.Alpha_f = -6.111 * 10**-6
        self.Alpha_p = 1.8e-6

        self.Fr = 0.975  # fission power factor
        self.Cpf = 1136  # fuel heat conductivity
        self.fuel_pitch = 1.26e-2
        self.Cpc = PropsSI("C", "P", Pressure, "Q", 0, "water")  # coolant conductivity
        # number of fuel rods
        # fuel diameter
        self.weight_of_a_fuel_rod = 2

        # volumes
        self.CoolantVcore = (
            (self.fuel_pitch**2 - np.pi * self.df**2 / 4)
            * self.active_height
            * self.N_rods
        )  # aikhane bvvhul ache
        self.area = np.pi * self.df * self.active_height * self.N_fue_rods

        self.FuelMass = (
            np.pi
            * (self.df - 2 * self.gas_gap) ** 2
            * self.active_height
            * self.N_fue_rods
            * self.rho_f
            / 4
        )
        self.Wc = 1889.68208

        """initial conditions"""
        self.NominalPower = Power
        self.Precursor = Precursor
        self.PowerRatio = 1
        self.ExternalReactivity = 0.0

        # control parameter will be used by control driving program

        if len(TemperatureFuel) != 5:
            raise ValueError("Here should be three initial conditions!")
        else:
            self.Tf1 = TemperatureFuel[0]
            self.Tf2 = TemperatureFuel[1]
            self.Tf3 = TemperatureFuel[2]
            self.Tf4 = TemperatureFuel[3]
            self.Tf5 = TemperatureFuel[4]

        if len(TemperatureModerator) != 10:
            raise ValueError("Here should be six initial conditions!")
        else:
            self.Tmo1 = TemperatureModerator[0]
            self.Tmo2 = TemperatureModerator[1]
            self.Tmo3 = TemperatureModerator[2]
            self.Tmo4 = TemperatureModerator[3]
            self.Tmo5 = TemperatureModerator[4]
            self.Tmo6 = TemperatureModerator[5]
            self.Tmo7 = TemperatureModerator[6]
            self.Tmo8 = TemperatureModerator[7]
            self.Tmo9 = TemperatureModerator[8]
            self.Tmo10 = TemperatureModerator[9]

        self.Tcl = Temp_coldleg
        self.Thl = Temp_hotleg
        self.Tup = Temp_upplerplenum
        self.Tlp = Temp_lowerplenum

        self.tempmod = np.sum(TemperatureModerator)
        self.tempfuel = np.sum(TemperatureFuel)
        self.piping_length = 37.00756804766576

        self.Pressure = Pressure

        """ necessary calculations """
        """mass"""
        self.core_center_factor = 0.05
        self.pipe_diameter = 850e-3
        self.RCP_position_factor = 0.5
        self.Mcl = (
            (self.piping_length - self.core_height)
            * self.RCP_position_factor
            * np.pi
            * self.pipe_diameter**2
            / 4
        ) * PropsSI("D", "P", self.Pressure, "Q", 0, "water")
        self.Mlp = (
            (
                (self.core_height - self.active_height)
                * self.d**2
                * PropsSI("D", "P", self.Pressure, "Q", 0, "water")
            )
            * self.core_center_factor
            / 4
        )
        self.Mup = (
            (
                (self.core_height - self.active_height)
                * self.d**2
                * PropsSI("D", "T", self.Tup, "Q", 0, "water")
            )
            * (1 - self.core_center_factor)
            / 4
        )
        self.Mhl = (
            (self.piping_length - self.core_height)
            * (1 - self.RCP_position_factor)
            * np.pi
            * self.pipe_diameter**2
            / 4
        ) * PropsSI("D", "P", self.Pressure, "Q", 0, "water")

        """this values still not known. """

        self.Mf = self.FuelMass / 5  # as there will  5  lumps
        self.Mmo = (
            self.CoolantVcore * self.CoolantDensity / 10
        )  # as there will 10  lumps

        # mass of the cold leg hot leg lower plenum
        """time const"""

        self.Lamda = self.total_delayed_const / (
            self.Beta1 / self.Lamda1
            + self.Beta2 / self.Lamda2
            + self.Beta3 / self.Lamda3
            + self.Beta4 / self.Lamda4
            + self.Beta5 / self.Lamda5
            + self.Beta6 / self.Lamda6
        )

        self.reactivity = 0

    """###########################################################################################"""
    """                                     Neutronics an Power                                   """
    """###########################################################################################"""

    def Reactivity(self, FuelTempSum: float, ModeratorTempSum):
        self.reactivity = (
            self.ExternalReactivity
            + self.Alpha_f * (self.tempfuel - FuelTempSum) / 5
            + self.Alpha_m * (self.tempmod - ModeratorTempSum) / 10
        )

    def DPowerRatio(self):
        dtdPP0 = (
            (self.reactivity - self.total_delayed_const) * self.PowerRatio
        ) / self.NGT + self.Lamda * self.Precursor
        return dtdPP0

    def DPrecoursor(self):
        dtdc = (
            self.total_delayed_const * self.PowerRatio / self.NGT
            - self.Lamda * self.Precursor
        )
        return dtdc

    """###########################################################################################"""
    """                                         Fuel Region                                       """
    """###########################################################################################"""

    def DTf1(self):
        dtdTf1 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.U * self.area * (self.Tmo1 - self.Tf1) / (self.Mf * self.Cpf)
        return dtdTf1

    def DTf2(self):
        dtdTf2 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.U * self.area * (self.Tmo3 - self.Tf2) / (self.Mf * self.Cpf)
        return dtdTf2

    def DTf3(self):
        dtdTf3 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.U * self.area * (self.Tmo5 - self.Tf3) / (self.Mf * self.Cpf)
        return dtdTf3

    def DTf4(self):
        dtdTf2 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.U * self.area * (self.Tmo7 - self.Tf4) / (self.Mf * self.Cpf)
        return dtdTf2

    def DTf5(self):
        dtdTf3 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.U * self.area * (self.Tmo9 - self.Tf5) / (self.Mf * self.Cpf)
        return dtdTf3

    """###########################################################################################"""
    """                                     MODERATOR REGION                                      """
    """###########################################################################################"""

    def DTmo1(self):

        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf1 - self.Tmo1) / (self.Mmo * self.Cpc)
            + (self.Tlp - self.Tmo1) / self.time_constmo
        )
        dtdTm01 = a + b
        return dtdTm01

    def DTmo2(self):

        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf1 - self.Tmo1) / (self.Mmo * self.Cpc)
            + (self.Tmo1 - self.Tmo2) / self.time_constmo
        )

        dtdTmo2 = a + b

        return dtdTmo2

    def DTmo3(self):

        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf2 - self.Tmo3) / (self.Mmo * self.Cpc)
            + (self.Tmo2 - self.Tmo3) / self.time_constmo
        )
        dtdTm03 = a + b
        return dtdTm03

    def DTmo4(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf2 - self.Tmo3) / (self.Mmo * self.Cpc)
            + (self.Tmo3 - self.Tmo4) / self.time_constmo
        )
        dtdTm04 = a + b
        return dtdTm04

    def DTmo5(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf3 - self.Tmo5) / (self.Mmo * self.Cpc)
            + (self.Tmo4 - self.Tmo5) / self.time_constmo
        )
        dtdTm05 = a + b
        return dtdTm05

    def DTmo6(self):

        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf3 - self.Tmo5) / (self.Mmo * self.Cpc)
            + (self.Tmo5 - self.Tmo6) / self.time_constmo
        )

        dtdTmo2 = a + b

        return dtdTmo2

    def DTmo7(self):

        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (self.Tmo6 - self.Tmo7) / self.time_constmo + self.U * self.area * (
            self.Tf4 - self.Tmo7
        ) / (self.Mmo * self.Cpc)
        dtdTm03 = a + b
        return dtdTm03

    def DTmo8(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf4 - self.Tmo7) / (self.Mmo * self.Cpc)
            + (self.Tmo7 - self.Tmo8) / self.time_constmo
        )
        dtdTm04 = a + b
        return dtdTm04

    def DTmo9(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf5 - self.Tmo9) / (self.Mmo * self.Cpc)
            + (self.Tmo8 - self.Tmo9) / self.time_constmo
        )
        dtdTm05 = a + b
        return dtdTm05

    def DTmo10(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        a = (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
        b = (
            self.U * self.area * (self.Tf5 - self.Tmo9) / (self.Mmo * self.Cpc)
            + (self.Tmo9 - self.Tmo10) / self.time_constmo
        )
        dtdTm06 = a + b
        return dtdTm06

    def DTcl(self, Temp_RCP: float):

        dtdTcl = self.Wc * (Temp_RCP - self.Tcl) / self.Mcl
        return dtdTcl

    def DThl(self):

        dtdTcl = self.Wc * (self.Tup - self.Thl) / self.Mhl
        return dtdTcl

    def DTup(self):
        self.Mup = (
            (
                (self.core_height - self.active_height)
                * self.d**2
                * PropsSI("D", "P", self.Pressure, "Q", 0, "water")
            )
            * (1 - self.core_center_factor)
            / 4
        )
        dtdTcl = self.Wc * (self.Tmo10 - self.Tup) / self.Mup
        return dtdTcl

    def DTlp(self):
        self.Mlp = (
            (
                (self.core_height - self.active_height)
                * self.d**2
                * PropsSI("D", "P", self.Pressure, "Q", 0, "water")
            )
            * self.core_center_factor
            / 4
        )
        dtdTcl = self.Wc * (self.Tcl - self.Tlp) / self.Mlp
        return dtdTcl

    def integrator(self, function, argsforfunction: list, intitial_cond, time_step):
        l = len(argsforfunction)

        if l == 0:
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


"""------------------------------------------------------------ Done---------------------------------------------------- """
