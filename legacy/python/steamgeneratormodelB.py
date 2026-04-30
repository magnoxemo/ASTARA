import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from Base import Base

class PrimaryFluidLump(Base):

    def __init__(
        self,
        Pressure: float,
        PrimaryLump1Temperature: float,
        PrimaryLump2Temperature: float,
    ):
        self.Mass = 1000  # mass in primary lump
        self.W = 18899.68208  # flow rate from the reactor
        self.Cp1 = PropsSI("C", "P", Pressure, "Q", 0, "water")
        self.Tp1 = PrimaryLump1Temperature
        self.Tp2 = PrimaryLump2Temperature

        self.number_of_utube = 3383
        self.inner_diameter = 0.019685
        self.outer_diameter = 0.022225
        self.length = 10.831712
        self.frist_lump_length = 1.05017116
        self.Upm = 12263.68  # Heat transfer co-efficient between steel and water

        self.Mass1 = (self.frist_lump_length / self.length) * self.Mass
        self.Mass2 = (1 - self.frist_lump_length / self.length) * self.Mass

        self.Spm1 = (
            np.pi
            * self.length
            * self.outer_diameter
            * self.number_of_utube
            * (self.frist_lump_length / self.length)
        )

        self.Spm2 = (
            np.pi
            * self.length
            * self.outer_diameter
            * self.number_of_utube
            * (1 - self.frist_lump_length / self.length)
        )
        self.timeconst = self.Mass / self.W

    def DTp1(self, InletPlenum: object, MetalLump: object):
        val1 = self.W * self.Cp1 * InletPlenum.Tpi
        val2 = self.Upm * self.Spm1 * (self.Tp1 - MetalLump.Tm1)
        val3 = self.W * self.Cp1 * self.Tp1
        dtdTp1 = (val1 - val2 - val3) / (self.Mass1 * self.Cp1)
        return dtdTp1

    def DTp2(self, MetalLump: object):
        val1 = self.W * self.Cp1 * self.Tp1
        val2 = self.Upm * self.Spm2 * (self.Tp2 - MetalLump.Tm2)
        val3 = self.W * self.Cp1 * self.Tp2
        dtdTp2 = (val1 - val2 - val3) / (self.Mass2 * self.Cp1)
        return dtdTp2



class MetalLump(Base):
    def __init__(self, MetalLump1Temperature: float, MetalLump2Temperature: float):

        self.Tm1 = MetalLump1Temperature
        self.Tm2 = MetalLump2Temperature

        self.Ums = (
            11186.216  # effective heat transfer co-efficient between water and steel
        )
        self.Cm = 460.547802  # specific heat of the metal
        self.number_of_utube = 3383
        self.inner_diameter = 0.019685
        self.outer_diameter = 0.022225
        self.metaldensity = 8050.095
        self.length = 10.831712
        self.frist_lump_length = 1.05017116

        self.Mass1 = (
            self.metaldensity
            * self.number_of_utube
            * self.length
            * np.pi
            * (self.outer_diameter**2 - self.inner_diameter**2)
            * (self.frist_lump_length / self.length)
            / 4
        )

        self.Mass2 = (
            self.metaldensity
            * self.number_of_utube
            * self.length
            * np.pi
            * (self.outer_diameter**2 - self.inner_diameter**2)
            * (1 - self.frist_lump_length / self.length)
            / 4
        )

        self.Sms1 = (
            np.pi
            * self.length
            * self.outer_diameter
            * self.number_of_utube
            * (self.frist_lump_length / self.length)
        )

        self.Sms2 = (
            np.pi
            * self.length
            * self.outer_diameter
            * self.number_of_utube
            * (1 - self.frist_lump_length / self.length)
        )

    def DTm1(self, PrimaryFluidLump: object, SecondaryFluidLump: object):

        val1 = (
            PrimaryFluidLump.Upm
            * PrimaryFluidLump.Spm1
            * (PrimaryFluidLump.Tp1 - self.Tm1)
        )
        val2 = self.Ums * self.Sms1 * (self.Tm1 - SecondaryFluidLump.Tsat)

        dtdTm1 = (val1 - val2) / (self.Mass * self.Cm)
        return dtdTm1

    def DTm2(self, SecondaryFluidLump: object):

        val1 = (
            PrimaryFluidLump.Upm
            * PrimaryFluidLump.Spm2
            * (PrimaryFluidLump.Tp2 - self.Tm2)
        )
        val2 = self.Ums * self.Sms2 * (self.Tm2 - SecondaryFluidLump.Tsat)

        dtdTm2 = (val1 - val2) / (self.Mass2 * self.Cm)
        return dtdTm2



class SecondaryFluidLump(Base):

    def __init__(
        self,
        inletTemperature: float,
        Pressure: float,
        FeedWaterFlowRate: float,
        EvaporationRate: float,
        SteamExitRate: float,
    ):

        self.Tfi = inletTemperature
        self.P = Pressure
        self.Tsat = PropsSI("T", "P", self.P, "Q", 1, "water")
        self.hfg = PropsSI("H", "P", self.P, "Q", 1, "water") - PropsSI(
            "H", "P", self.P, "Q", 0, "water"
        )
        self.Vfg = PropsSI("V", "P", self.P, "Q", 0, "water") - PropsSI(
            "V", "P", self.P, "Q", 1, "water"
        )
        self.Msw = 1000  # mass of water in the secondary lump
        self.Mss = 1002  # mass of the Saturated steam in the secondary lump

        self.Wfi = FeedWaterFlowRate
        self.Wsg = EvaporationRate
        self.Wso = SteamExitRate
        from scipy.interpolate import interp1d

        """ Gradient Function Calculation """

        pressure = np.linspace(1e3, 2.7e7, num=20000)
        Hf = []
        Hg = []

        for i in pressure:
            Hf.append(PropsSI("H", "P", i, "Q", 0, "water"))
            Hg.append(PropsSI("H", "P", i, "Q", 1, "water"))

        TempGrad = np.gradient(pressure)
        HfGrad = np.gradient(Hf)
        HgGrad = np.gradient(Hg)

        k1 = HfGrad / TempGrad
        k2 = HgGrad / TempGrad

        self.dHfdTstat = interp1d(pressure, k1)
        self.dHgdTstat = interp1d(pressure, k2)

    def DPressure(self, MetalLump: object, PrimaryFluidLump: object):

        val1 = (
            MetalLump.Ums
            * MetalLump.Sms1
            * (MetalLump.Tm1 + MetalLump.Tm2 - self.Tsat * 2)
        )
        val2 = self.Wfi * PrimaryFluidLump.Cp2 * self.Tfi
        val3 = PropsSI("H", "P", self.P, "Q", 1, "water") * (self.Wsg - self.Wso)
        val4 = PropsSI("H", "P", self.P, "Q", 0, "water") * (
            self.Wfi - self.Wsg
        ) + self.Wso * PropsSI("H", "P", self.P, "Q", 1, "water")

        dtdP = (val1 + val2 - val3 - val4) / (
            self.Msw * self.dHfdTstat(self.P) + self.Mss * self.dHgdTstat(self.P)
        )

        return dtdP

    def DMsw(self):

        dtdMsw = self.Wfi - self.Wsg

        return dtdMsw

    def DMss(self):

        dtdMss = self.Wsg - self.Wso

        return dtdMss


