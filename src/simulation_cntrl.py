import keyboard 
import sys 
import thread
import os 

from UTSG_modeling import UTSG
from pump import pump
from Reactor_model import Reactor

print("Press Q key from the keyboard to stop the simulation ")

while True:
  if keyboard.is_pressed('q') or keyboard.is_pressed('Q'):
    break
  else:
     #rest of the code +multi threading will be put here
     pass 
