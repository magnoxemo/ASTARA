def step_integrator(function,initial_condition=(x0,y0),stepsize):

"""the method we will be using here is runge kutta 4th order method
this is a step wise function.so initial condition for nth step will
be Y(n-1) value.

function           --f(x,y) -> the differential function to be integrated shape y'=f(x,y)
initial condition  --(x0,y0)-> this will be a tuple x:independent variable y:dependent variable  
step size          --h      -> a float value """
	stepsize=h
	
	b1=function(x0,y0)
	b2=function(x0+h*b1/2,y0+h*b1/2)
	b3=function(x0+h*b2/2,y0+h*b2/2)
	b4=function(x0+h*b3,y0+h*b1)
	
	return y0+h*(b1+2*b2+2*b3+b4)/6
