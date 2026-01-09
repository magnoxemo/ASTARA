import numpy as np

# --------------------------Reactor constants------------------------
beta = np.array(
    [0.000209, 0.001414, 0.001309, 0.002727, 0.000925, 0.000314]
)  # neutron group constants
BETAt = np.sum(beta)  # total neutron group constants

LAmda = np.array(
    [0.0125, 0.0308, 0.114, 0.307, 1.19, 3.19]
)  # neutron group decay constants

CO = 4625.71
Pngt = 17.9 * 10**-6  # prompt neutron generation time

ALPHAf = -1.1 * 10**-5  # Fuel Coefficient of Reactivity
ALPHAc = -2 * 10**-4  # Moderator Coefficient of Reactivity
ALPHAp = -1 * 10**-6

Mass_fuel = 222739  # mass of fuel
Cpf = 0.059
FR = 0.974  # fraction of total power generated in the fuel
Vup = 1376  # Volume of the coolant in the upper pennulum
Vlp = 1791  # volume of the coolant in the lower pennulum
Vhl = 1000  # Volume of the coolant in the hot leg pipings
Vcl = 2000  # Volume of the coolant in the cold leg pipings
Vmo = 540  # Volume of the coolant in the core
DENSrno = 45.71
Cprno = 1.39

Bcl = -0.0631
Bhl = -0.0741

Tpo0 = 539.5
Tf10 = 1488.75
Trno10 = 548.303
Trno20 = 557.106
Tf20 = 1506.38
Trno30 = 565.909
Trno40 = 574.711
Tf30 = 1523.96
Trno50 = 583.514
Trno60 = 592.5
Thl0 = 592.5
Tup0 = 592.5
Tlp0 = 539.5
Tc10 = 539.5

ROHexO = 0
POWERi = 691244.4199
Tset0 = 566
TOUset = 30
Taves0 = 566
TOUla1 = 10
TOUla2 = 5
TOUle = 80
AUXCOO = 0
Tclp0 = 539.5
Thlp0 = 592.5
TOUrtd = 4
React = 0.003375
Pcor0 = 1
Pcor10 = 3436 * 10**9
Wpth = 1.57559e8

"""End of the reactor costants """

"""constant calculation of the core """

Pcor1 = (Pcor10 * 3.41) / (3600 * 3)
Hfm = 200 / 3600
Afm = 59900 / 3
Wpt = Wpth / 3600

s = 0

for i in range(len(beta)):
    s += beta[i] / LAmda[i]
LANDA = BETAt / s

Mf = Mass_fuel / 3
Mmo = Vmo * DENSrno
Mup = Vup * DENSrno
Mlp = Vlp * DENSrno
Mhl = Vhl * DENSrno
Mcl = Vcl * DENSrno

TOUmo = Mmo / (2 * Wpt)
TOUmo = Mup / Wpt
TOU1p = Mlp / Wpt
TOUhl = Mhl / Wpt
TOUcl = Mhl / Wpt

"""-----------------------------initail conditions-------------------------"""
###Fuel temp in 3 nodes needs to be determined
Tf1 = 3
Tf2 = 3
Tf3 = 34
Tmo1 = 1
Tmo2 = 232
Tmo3 = 2
Tmo4 = 1
Tmo5 = 2
Tmo6 = 3


Tlp = 121
C = 121
Pcor = 2

"""code which will run on a loop basically formula of differential parts and loop integration  """

ROH = (
    ROHexO
    + ALPHAf * (Tf1 + Tf2 + Tf3 - Tf10 - Tf20 - Tf30) / 3
    + ALPHAc
    * (
        Tmo1
        + Tmo2
        + Tmo3
        + Tmo4
        + Tmo5
        + Tmo6
        - Trno10
        - Trno20
        - Trno30
        - Trno40
        - Trno50
        - Trno60
    )
    / 6
)
dPcor_dt = (ROH - BETAt) * Pcor / Pngt + LANDA * C
dC_dt = BETAt * Pcor / Pngt - LANDA * C


dTf1_dt = (FR * Pcor1 * Pcor) / (Mass_fuel * Cpf) + (Hfm * Afm) * (Tmo1 - Tf1) / (
    Mass_fuel * Cpf
)
dTmo1_dt = (
    ((1 - FR) * Pcor1 * Pcor) / (Mmo * Cprno)
    + (Hfm * Afm) * (Tf1 - Tmo1) / (Mmo * Cprno)
    + (Tlp - Tmo1) / TOUmo
)
dTmo2_dt = (
    ((1 - FR) * Pcor1 * Pcor) / (Mmo * Cprno)
    + (Hfm * Afm) * (Tf1 - Tmo1) / (Mmo * Cprno)
    + (Tmo1 - Tmo2) / TOUmo
)

dTf2_dt = (FR * Pcor1 * Pcor) / (Mass_fuel * Cpf) + (Hfm * Afm) * (Tmo3 - Tf2) / (
    Mass_fuel * Cpf
)
dTmo3_dt = (
    ((1 - FR) * Pcor1 * Pcor) / (Mmo * Cprno)
    + (Hfm * Afm) * (Tf2 - Tmo3) / (Mmo * Cprno)
    + (Tmo2 - Tmo3) / TOUmo
)
dTmo4_dt = (
    ((1 - FR) * Pcor1 * Pcor) / (Mmo * Cprno)
    + (Hfm * Afm) * (Tf2 - Tmo3) / (Mmo * Cprno)
    + (Tmo3 - Tmo4) / TOUmo
)

dTf2_dt = (FR * Pcor1 * Pcor) / (Mass_fuel * Cpf) + (Hfm * Afm) * (Tmo5 - Tf3) / (
    Mass_fuel * Cpf
)
dTmo5_dt = (
    ((1 - FR) * Pcor1 * Pcor) / (Mmo * Cprno)
    + (Hfm * Afm) * (Tf3 - Tmo5) / (Mmo * Cprno)
    + (Tmo4 - Tmo5) / TOUmo
)
dTmo6_dt = (
    ((1 - FR) * Pcor1 * Pcor) / (Mmo * Cprno)
    + (Hfm * Afm) * (Tf3 - Tmo5) / (Mmo * Cprno)
    + (Tmo5 - Tmo6) / TOUmo
)
