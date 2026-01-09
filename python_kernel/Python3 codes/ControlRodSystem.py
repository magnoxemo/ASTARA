class ControlRodSystem():

    def __init__(self,time:float):

      '''
      Insertation speed and the withdrawal speed can be adjusted if wishes to other wise the operation will be by default.
      the StepReactivity value depends on the type of reactor core.So it needs to be adjusted if someone not wishes to VVER-1200
      for their simulation.
      '''
        self.Insertationspeed=0.5
        self.Withdrawalspeed=0.4
        self.StepReactivity=0.002
        self.time=time

    def DRouex(self,Mode:str):
        if Mode=='insert':
            dtdRouex=self.Insertationspeed*self.StepReactivity
            return dtdRouex
        elif Mode=='withdraw':
            dtdRouex=self.Withdrawalspeed*self.StepReactivity
            return dtdRouex
        elif Mode=='Hover':
            dtdRouex=0
            return dtdRouex
        else:
            message=" Control rod operation mode can't be "+ Mode
            raise AttributeError(message)
