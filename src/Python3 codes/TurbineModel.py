from CoolProp.CoolProp import PropsSI
import numpy as np
from Base import Base

class Turbine:

    def __init__(
        self,
        Throttle_valve: object,
        MoistureSeperator: object,
        Reheater: object,
        HighPressureTurbine: object,
        LowPressureTurbine: object,
    ):
        self.ThrottleValve = Throttle_valve
        self.MoistureSeperator = MoistureSeperator
        self.Reheater = Reheater
        self.HighPressureTurbine = HighPressureTurbine
        self.LowPressureTurbine = LowPressureTurbine


class ThrottleValve:
    def __init__(
        self,
        Area_main: float,
        Area_secondary: float,
        co_efficient_main: float,
        co_efficient_secondary: float,
        pos_main_valve: float,
        pos_second_valve: float,
        enthalpy_at_main_throttle_valve: float,
    ):
        """Desgin Parameters"""
        self.A_m = Area_main
        self.A_s = Area_secondary
        self.Cf_m = co_efficient_main
        self.Cf_s = co_efficient_secondary

        self.W_utsg = 4 * 1329.6
        self.SteamTemp = 265  # temperature will be in Kelvin
        self.SteamPressure = 5e6
        # pressure will be in Pa
        """ SteamTemp and the SteamPressure will be from the steam generator  
            So once the constructor runs  in the integrator loop I need to update the 
            fluid properties 
            
            Things to keep in mind 
            
            1.Temp will be in Kelvin 
            2.The pressure will be in Pascal

        """

        self.h_s = PropsSI("H", "T", self.SteamTemp, "P", self.SteamPressure, "water")

        """ control elements """
        if (pos_main_valve > 100 and pos_main_valve < 0) or (
            pos_second_valve > 100 and pos_second_valve < 0
        ):
            raise ValueError(
                "opening percentage of the  Main and secondarythrottle valve must stay in between 0 and 100"
            )
        else:
            self.Pos_main_v = pos_main_valve
            self.pos_second_v = pos_second_valve

        self.Wm = self.Cf_m * self.A_m * self.W_utsg * self.Pos_main_v
        self.W2nd = self.Cf_s * self.A_s * self.W_utsg * self.pos_second_v

        # position of the main and second valve means how much the percentage of the valve is opened
        # Here is Pos_main_V=100 then it's fully open

    def Wmain(self):
        self.Wm = self.Cf_m * self.A_m * self.W_utsg * self.Pos_main_v
        self.W2nd = self.Cf_s * self.A_s * self.W_utsg * self.pos_second_v


class NozzleChest(Base):

    def __init__(
        self,
        Effective_volume_of_nozzle_Chest: float,
        steam_pressure_chest: float,
        TempatNozzleChest: float,
        Kc_hp: float,
        Callender_const1: float,
        Callender_const2: float,
    ):
        """k1 and k2 are constants in the Callender’s emperical
        equation relating pressure, density and enthalpy of superheated steam

        Kchp is a constant determined at the initial period of the powerplant
        """
        self.k1 = Callender_const1
        self.k2 = Callender_const2
        self.Kchp = Kc_hp

        """Data loading from the CoolProP"""

        self.SteamTemp = 265  # temperature will be in Kelvin
        self.SteamPressure = 5e6
        self.h_sd = PropsSI("H", "T", self.SteamTemp, "P", self.SteamPressure, "water")
        # enthalpy of the steam drum

        self.Pc = steam_pressure_chest
        self.SteamTempatChest = TempatNozzleChest
        self.Vc = Effective_volume_of_nozzle_Chest
        self.hc = PropsSI("H", "T", self.SteamTempatChest, "P", self.Pc, "water")
        self.rou_c = PropsSI("D", "T", self.SteamTempatChest, "P", self.Pc, "water")

    def Whp(self, Reheater: object, HighPressureTurbine: object):
        self.Whp1 = self.Kchp * np.sqrt(
            self.Pc * self.rou_c - Reheater.Pressue * HighPressureTurbine.rou_exit
        )

    def P_c(self):
        self.Pc = PropsSI("P", "T", self.SteamTempatChest, "Q", 0, "water")

    def Dh_c(self, ThrottleValve: object):
        self.Wmain = ThrottleValve.Wm
        dtdhc = (
            ((self.Wmain * self.h_sd - self.Whp1 * self.hc) / (self.rou_c * self.Vc))
            + (self.Pc / self.rou_c**2) * self.Drou_c()
        ) / (1 - self.k1)

        return dtdhc

    def Drou_c(self, ThrottleValve: object):
        self.Wmain = ThrottleValve.Wm
        dtdrou_c = (self.Wmain - self.Whp1) / self.Vc
        return dtdrou_c


    """ As we can get the transient feed back of the enthalpy and pressure and density 
    we can find out the temperature of the steam using cool prop properties """


class MoistureSeperator:
    def __init__(
        self,
        enthalpy_hpex: float,
        flow_rate_of_the_steam_from_the_steamseperator: float,
        Temperature: float,
    ):

        self.Whpex = flow_rate_of_the_steam_from_the_steamseperator  # 229.2 kg/s
        self.Temp = Temperature
        self.h_hpex = enthalpy_hpex
        self.hf = PropsSI("H", "T", self.Temp, "Q", 0, "water")
        self.hg = PropsSI("H", "T", self.Temp, "Q", 1, "water")
        self.hfg = self.hg - self.hf

        # variable
        self.W_mss = (self.h_hpex - self.hf) * self.Whpex / (self.hfg)

    def Wmsw(self, HighPressureTurbine: object, NozzleChest: object):
        self.Wbhp = HighPressureTurbine.Wbhp
        self.Whp1 = NozzleChest.Whp()
        self.W_msw = (self.W_hp1 - self.Wbhp) - (
            (self.h_hpex - self.hf) * self.Whpex / (self.hfg)
        )

    def Wmss(self):
        self.W_mss = (
            (self.h_hpex - self.hf) * self.Whpex / (self.hfg)
        )  # that will go to the reheater


class Reheater(Base):

    def __init__(
        self,
        MoistureSpererator: object,
        VolumeOftheReheater: float,
        time_const_flowrate: float,
        time_const_heating: float,
        flow_rate_to_2nd_heater: float,
        Heater_temp: float,
        heat_transfer_coefficient: float,
        Gas_const: float,
        Throttlechestflowrate: float,
    ):
        # self.W_mss=MoistureSeperator.W_mss

        self.hr = 2920.2e3  # enthalpy at the steam reheater exit
        self.Vr = VolumeOftheReheater
        self.Density = 2.4
        self.Heat = 4.365e7  # reheater heat
        self.Temperature = PropsSI("T", "D", self.Density, "Q", 1, "water")
        self.hmss = PropsSI("H", "T", self.Temperature, "Q", 0, "water")
        self.Pressure = PropsSI("P", "T", self.Temperature, "Q", 0, "water")
        self.Kclp = 0.12
        self.Wlp = self.Kclp * np.sqrt(self.Pressure * self.Density)
        self.k1 = 0.01

        self.w_2nd = Throttlechestflowrate
        self.W_ro = flow_rate_to_2nd_heater  # for initial_cond

        self.tau_1 = time_const_flowrate
        self.tau_2 = time_const_heating
        self.T_steam = self.Temperature
        self.T_r = Heater_temp

        self.H = heat_transfer_coefficient
        self.R = Gas_const

    def DDensity(self):
        dtdDensity = (MoistureSeperator.W_mss - self.Wlp) / self.Vr
        return dtdDensity

    def Dhr(self):
        dtdhr = (
            (
                (self.Heat + MoistureSeperator.W_mss * self.hmss - self.Wlp * self.hr)
                / (self.Density * self.Vr)
            )
            + (
                self.Pressure
                * (MoistureSeperator.W_mss - self.Wlp)
                / (self.Density**2 * self.Vr)
            )
        ) * (1 / (1 - self.k1))

        return dtdhr

    def DWro(self):

        dtdWro = (self.w_2nd - self.W_ro) / self.tau_1

        return dtdWro

    def DQr(self):

        dtdQr = (
            (self.T_steam - self.T_r) * (self.w_2nd + self.W_ro) * self.H - 2 * self.Q
        ) / (2 * self.tau_2)
        self.T_r = self.P / (self.R * self.rou_r)

        return dtdQr



class HighPressureTurbine(Base):
    def __init__(
        self,
        inlet_flow_rate: float,
        PressureIn: float,
        PressureOut: float,
        time_const: float,
        HP_co_efficient: float,
    ):

        self.Whpin = inlet_flow_rate
        self.PressureIN = PressureIn
        self.PressureOUT = PressureOut
        self.Whpex = (1 - self.C) * self.Whpin
        self.WhpexMax = 299.2
        self.Wbhp = self.C * self.Whpin
        self.Tau = time_const
        self.Ws = self.Wbhp
        self.EfficiencyCorrectionFactor = 0.01
        self.Efficiency = 0.86

        self.Hc = PropsSI("H", "P", self.PressureIN, "Q", 1, "water")
        self.Hex = PropsSI("H", "P", self.PressureOUT, "Q", 0.9, "water")

        self.C = HP_co_efficient

    def DWhpex(self):

        dtdWhpex = ((self.Whpin - self.Wbhp) - self.Whpex) / self.Tau
        self.Wbhp = self.C * self.Whpin

        return dtdWhpex

    def StageEfficiency(self):
        self.neu = (
            self.Efficiency
            * ((self.Whpex / self.WhpexMax) - self.EfficiencyCorrectionFactor)
            / ((self.Whpex / self.WhpexMax) + 1 - self.EfficiencyCorrectionFactor)
        )

    def Torque(self):

        self.Torqu = self.neu * self.Whpin * (self.Hc - self.Hex) / (2 * np.pi**2 * 120)
        # as the angular frequency is 120pi



class LowPressureTurbine(Base):
    def __init__(
        self,
        inlet_flow_rate: float,
        exit_flow_rate_to_MS: float,
        PressureIN: float,
        PressureOUT: float,
        time_const: float,
        LP_co_efficient: float,
    ):

        self.Wlpex = exit_flow_rate_to_MS  # goes to the condenser
        self.WlpexMax = 238.8
        self.Wlpin = inlet_flow_rate
        self.PressureIN = PressureIN
        self.PressureOUT = PressureOUT
        self.C = LP_co_efficient  # goes to the heater
        self.Wblp = self.C * self.Wlpin
        self.Tau = time_const  # this comes from the moisture seperator and reheater
        self.C = LP_co_efficient
        self.EfficiencyCorrectionFactor = 0.01
        self.Efficiency = 0.86
        self.omega = 120 * np.pi

        self.Hc = PropsSI("H", "P", self.PressureIN, "Q", 1, "water")
        self.Hex = PropsSI("H", "P", self.PressureOUT, "Q", 0.9, "water")

    def DWhpex(self):

        dtdwlpex = ((self.Wlpin - self.Wblp) - self.Wlpex) / self.Tau
        self.Wblp = self.C * self.Wlpin

        return dtdwlpex

    def StageEfficiency(self):
        self.neu = (
            self.Efficiency
            * ((self.Wlpex / self.WlpexMax) - self.EfficiencyCorrectionFactor)
            / ((self.Wlpex / self.WlpexMax) + 1 - self.EfficiencyCorrectionFactor)
        )

    def Torque(self):

        self.Torqu = (
            self.neu * self.Wlpin * (self.Hc - self.Hex) / (2 * np.pi * self.omega)
        )
        # as the angular frequency is 120pi


class LowPressureHeater(Base):
    def __init__(
        self,
        FeedWaterFlowRate: float,
        ExitFlowRateHP: float,
        ExitFlowRateMS: float,
        ExitFlowRateRH: float,
        timeConst1: float,
        timeconst2: float,
        PressureCondenser: float,
        TemperatureCondenser: float,
    ):

        self.Hco = PropsSI("H", "P", PressureCondenser, "T", TemperatureCondenser)
        self.Hfw = 425.4
        self.Wfw = FeedWaterFlowRate
        self.Wblp = ExitFlowRateHP
        self.Wmsw = ExitFlowRateMS
        self.Wro = ExitFlowRateRH
        self.timeconst1 = timeConst1
        self.timeconst2 = timeconst2
        self.H = 1.339e6

    def Dhfw(self):
        dtdhfw = (
            self.H * (self.Wblp + self.Wmsw + self.ro) / (self.timeconst1 * self.Wfw)
            + (self.Hco - self.Hfw) / self.timeconst1
        )
        return dtdhfw

    def DWhpd(self):
        dtdWhpd = self.Whp - self.Wmsw + self.Wro + self.Whpd

        return dtdWhpd




class HighPressureHeater(Base):
    def __init__(
        self,
        FeedWaterFlowRate: float,
        Whpd: float,
        ExitFlowRatelp: float,
        ExitFlowRateRH: float,
        timeConst: float,
        PressureCondenser: float,
        TemperatureCondenser: float,
    ):

        self.Hco = PropsSI("H", "P", PressureCondenser, "T", TemperatureCondenser)
        self.Hfw = 425.4
        self.Wfw = FeedWaterFlowRate
        self.Wblp = ExitFlowRatelp
        self.Whpd = Whpd
        self.timeconst = timeConst
        self.H = 1.339e6

    def Dhfw(self):
        dtdhfw = (
            self.H * (self.Whpd + self.Wblp) / (self.timeconst * self.Wfw)
            + (self.Hco - self.Hfw) / self.timeconst
        )
        return dtdhfw



# modeling done
