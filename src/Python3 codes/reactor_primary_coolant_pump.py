form Base import Base
class reactor_primary_coolant_pump(Base):

    def __init__(
        self,
        T_po: float,
        pipe_diameter: float,
        piping_length: float,
        friction_factor: float,
        moment_of_innertia: float,
        pump_speed: float,
        flow_rate: float,
        power_delivered_to_shaft: float,
        const: list,
    ):

        self.A = np.pi * pipe_diameter**2 / 4
        self.k = friction_factor
        self.L = piping_length
        self.g = 9.8
        self.water_density = PropsSI("D", "T", T_po, "Q", 0, "water")
        self.const = const
        self.N = pump_speed
        self.Q = flow_rate
        self.Pd = power_delivered_to_shaft
        self.I = moment_of_innertia
        self.Head_update()

    def Head_update(self):
        sum = 0
        for i in range(len(self.const)):
            sum = sum + self.const[i] * self.Q**i
        self.H = sum

    def DQ(self):

        Dq = self.A * self.g * (self.H - self.k * self.Q) / self.L

        return Dq

    def DN(self, T_po: float):
        self.water_density = PropsSI("D", "T", T_po, "Q", 0, "water")
        DNp = (self.Pd - self.g * self.water_density * self.Q * self.H) / (
            self.N * self.I * 4 * np.pi**2
        )

        return DNp

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
