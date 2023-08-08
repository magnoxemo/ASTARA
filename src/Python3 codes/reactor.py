import numpy as np 
import matplotlib.pyplot as plt 
import scipy as sp
import threading as th

""" Reactor modeling """

class CoreInside():
    def __init__(self):
        self.diameter=36.2712
        self.core_height=43.8912

        #group const
        self.Beta1=0.000209
        self.Beta2=0.001414
        self.Beta3=0.001414
        self.Beta4=0.00272 
        self.Beta5=0.00092
        self.Beta6=0.00689

        self.total_delayed_const=self.Beta1+self.Beta2+self.Beta3+self.Beta4\
            +self.Beta5+self.Beta6
        
        #decay const
        self.Lamda1=0.0125
        self.Lamda2=0.0308
        self.Lamda3=0.114
        self.Lamda4=0.307
        self.Lamda5=1.19
        self.Lamda6=3.19

        self.NGT=17.9*10**-6
        self.Alpha_m=-1.1192*10**-7
        self.Alpha_f=6.111*10**-6

        self.Fr=0.974       #fission power factor 
        self.Cpf=247.5052   #coolant conductivity 
        self.Cpmo=5819.652  #fuel heat conductivity 
        #volumes 

        self.Vup=38.9813
        self.Vlp=50.7091
        self.Vhl=28.3168
        self.Vcl=56.63369
