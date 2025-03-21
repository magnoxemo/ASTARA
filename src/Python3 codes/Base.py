class Base:
    
    def integrator(self, function, argsforfunction: None, intitial_cond, time_step):
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
