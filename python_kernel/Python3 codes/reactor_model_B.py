"""########################################################################################################"""

"""                                      Less Lump Model                                                   """
"""########################################################################################################"""

import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from tqdm import tqdm
form Base import Base

class CoreInside(Base):
    def __init__(
        self,
        TemperatureFuel: list,
        TemperatureModerator: list,
        Power: float,
        Precursor: float,
        Pressure: float,
    ):

        self.diameter = 3.62712
        self.core_height = 4.38912
        self.CoolantDensity = PropsSI("D", "P", Pressure, "Q", 0, "water")
        self.area = 14855.196091203 / 3
        self.h = 1134  # needs more exect value conductivity

        # group const
        self.Beta1 = 0.000243
        self.Beta2 = 0.001363
        self.Beta3 = 0.001203
        self.Beta4 = 0.002605
        self.Beta5 = 0.000829
        self.Beta6 = 0.000166

        self.total_delayed_const = (
            self.Beta1 + self.Beta2 + self.Beta3 + self.Beta4 + self.Beta5 + self.Beta6
        ) * 1.1

        # decay const
        self.Lamda1 = 0.0127
        self.Lamda2 = 0.0317
        self.Lamda3 = 0.115
        self.Lamda4 = 0.311
        self.Lamda5 = 1.40
        self.Lamda6 = 3.87

        self.NGT = 2 * 10**-5
        self.Alpha_m = -1.1192 * 10**-7
        self.Alpha_f = -6.111 * 10**-6

        self.Fr = 0.974  # fission power factor
        self.Cpc = PropsSI("C", "P", Pressure, "Q", 0, "water")  # coolant conductivity
        self.Cpf = 1136  # fuel heat conductivity
        # volumes
        self.CoolantVcore = 15.2911

        self.FuelMass = 101032.711
        self.Wc = 18899.68208

        """initial conditions"""
        self.NominalPower = Power
        self.Precursor = Precursor
        self.PowerRatio = 1
        self.ExternalReactivity = 0.00000

        # control parameter will be used by control driving program

        if len(TemperatureFuel) != 3:
            raise ValueError("Here should be three initial conditions!")
        else:
            self.Tf1 = TemperatureFuel[0]
            self.Tf2 = TemperatureFuel[1]
            self.Tf3 = TemperatureFuel[2]

        if len(TemperatureModerator) != 6:
            raise ValueError("Here should be six initial conditions!")
        else:
            self.Tmo1 = TemperatureModerator[0]
            self.Tmo2 = TemperatureModerator[1]
            self.Tmo3 = TemperatureModerator[2]
            self.Tmo4 = TemperatureModerator[3]
            self.Tmo5 = TemperatureModerator[4]
            self.Tmo6 = TemperatureModerator[5]

        self.tempmod = np.sum(TemperatureModerator)
        self.tempfuel = np.sum(TemperatureFuel)
        """ necessary calculations """
        """mass"""
        self.Mf = self.FuelMass / 3  # as there will three lumps
        self.Mmo = self.CoolantVcore * self.CoolantDensity / 3

        """time const"""
        self.time_constmo = self.Mmo / (self.Wc * 2)

        self.Lamda = self.total_delayed_const / (
            self.Beta1 / self.Lamda1
            + self.Beta2 / self.Lamda2
            + self.Beta3 / self.Lamda3
            + self.Beta4 / self.Lamda4
            + self.Beta5 / self.Lamda5
            + self.Beta6 / self.Lamda6
        )

        self.reactivity = 0

    def Reacivity(self, FuelTempSum: float, ModeratorTempSum):

        self.reactivity = (
            self.ExternalReactivity
            + self.Alpha_f * (self.tempfuel - FuelTempSum) / 3
            + self.Alpha_m * (self.tempmod - ModeratorTempSum) / 6
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

    def DTf1(self):
        dtdTf1 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.h * self.area * (self.Tmo1 - self.Tf1) / (self.Mf * self.Cpf)
        return dtdTf1

    def DTf2(self):
        dtdTf2 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.h * self.area * (self.Tmo3 - self.Tf2) / (self.Mf * self.Cpf)
        return dtdTf2

    def DTf3(self):
        dtdTf3 = self.Fr * self.NominalPower * self.PowerRatio / (
            self.Mf * self.Cpf
        ) + self.h * self.area * (self.Tmo5 - self.Tf3) / (self.Mf * self.Cpf)
        return dtdTf3

    def DTm01(self, LowerPlenum: object):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        dtdTm01 = (
            (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
            + self.h * self.area * (self.Tf1 - self.Tmo1) / (self.Mmo * self.Cpc)
            + (LowerPlenum.Tlp - self.Tmo1) / self.time_constmo
        )
        return dtdTm01

    def DTm02(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        dtdTm02 = (
            (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
            + self.h * self.area * (self.Tf1 - self.Tmo1) / (self.Mmo * self.Cpc)
            + (self.Tmo1 - self.Tmo2) / self.time_constmo
        )

        return dtdTm02

    def DTm03(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        dtdTm03 = (
            (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
            + self.h * self.area * (self.Tf2 - self.Tmo3) / (self.Mmo * self.Cpc)
            + (self.Tmo2 - self.Tmo3) / self.time_constmo
        )

        return dtdTm03

    def DTm04(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        dtdTm04 = (
            (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
            + self.h * self.area * (self.Tf2 - self.Tmo3) / (self.Mmo * self.Cpc)
            + (self.Tmo3 - self.Tmo4) / self.time_constmo
        )
        return dtdTm04

    def DTm05(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        dtdTm05 = (
            (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
            + self.h * self.area * (self.Tf3 - self.Tmo5) / (self.Mmo * self.Cpc)
            + (self.Tmo4 - self.Tmo5) / self.time_constmo
        )
        return dtdTm05

    def DTm06(self):
        self.time_constmo = self.Mmo / (self.Wc * 2)
        dtdTm06 = (
            (1 - self.Fr) * self.NominalPower * self.PowerRatio / (self.Mmo * self.Cpc)
            + self.h * self.area * (self.Tf3 - self.Tmo5) / (self.Mmo * self.Cpc)
            + (self.Tmo5 - self.Tmo6) / self.time_constmo
        )
        return dtdTm06

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


class Hotleg:
    def __init__(self, TempHotLeg: float):
        self.Thl = TempHotLeg
        self.CoolantVhl = 28.3168
        self.Wc = 18899.68208
        self.CoolantDensity = 732.18278
        self.Mhl = self.CoolantVhl * self.CoolantDensity
        self.time_consthl = self.Mhl / self.Wc

    def DThl(self, UpperPlenum: object):

        dtdThl = (UpperPlenum.Tup - self.Thl) / self.time_consthl
        return dtdThl

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


class ColdLeg:
    def __init__(self, TempColdLeg: float) -> None:
        self.Tcl = TempColdLeg
        self.Wc = 18899.68208
        self.CoolantDensity = 732.18278
        self.CoolantVcl = 56.63369
        self.Mcl = self.CoolantVcl * self.CoolantDensity
        self.time_constcl = self.Mcl / self.Wc

    def DTcl(self, TempRCPexit: float):
        """temperature exit will accessed through objects. but right now it's initilized as
        a float data type"""
        dtdTcl = (TempRCPexit - self.Tcl) / self.time_constcl
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


class UpperPlenum:
    def __init__(self, TempUpperPlenum: float):
        self.Tup = TempUpperPlenum
        self.Wc = 18899.68208
        self.CoolantDensity = 732.18278
        self.CoolantVup = 38.9813
        self.Mup = self.CoolantDensity * self.CoolantVup
        self.time_constup = self.Mup / self.Wc

    def DTup(self, ReactorCoreInside: object):
        dtdTup = (ReactorCoreInside.Tmo6 - self.Tup) / self.time_constup
        return dtdTup



class LowerPlenum(Base):
    def __init__(self, TempLowerPlenum: float):
        self.Tlp = TempLowerPlenum
        self.Wc = 18899.68208
        self.CoolantDensity = 732.18278
        self.CoolantVlp = 50.7091
        self.Mlp = self.CoolantVlp * self.CoolantDensity
        self.time_constlp = self.Mlp / self.Wc

    def DTlp(self, ColdLeg: object):
        dtdTlp = (ColdLeg.Tcl - self.Tlp) / self.time_constlp
        return dtdTlp

