# This function allows to use ranges, just as the built-in function range(), but with float arguments.
# Taken from "http://code.activestate.com/recipes/66472-frange-a-range-function-with-float-increments/"
def frange(start, end=None, inc=None):
	if end == None:
		end = start + 0.0
		start = 0.0

	if inc == None:
		inc = 1.0

	L = []
	while 1:
		next = start + len(L) * inc
		if inc > 0 and next >= end:
			break
		elif inc < 0 and next <= end:
			break
		L.append(next)
		
	return L


# The function "func" is the sigmoid fit of the density data points.
def func(r, ro, R, d):
	return (ro/2)*(1-np.tanh(2*(r-R)/d))
			
# Searching for initial guess for the first parameter of the fit
def guess(y):
	m = max(y)
	m2 = min(y)
	third = (m - m2)/3
	bottom = m - third
	tot = 0
	count = 0
	for i in range(0, len(y)):
		if y[i] >= bottom:
			tot = tot + y[i]
			count = count + 1
	At = tot / count
	return At


# FUNCTION TO DETERMINE BULK VALUE OF WATER DENSITY
# It takes 3 arguments: 
# input: all the data of the density map of one system
# zcoord: the z-coordinates of the slices of the system
# func: sigmoidal function used to fit the density data of one slice
def waterbulk(input, zcoord, func):
	x = input[0]
	densmean = np.zeros(200)
	i = 0
	start = 0
	for slice in range(1,200):
		ymean = input[slice]
		At = guess(ymean)					
		popt, pcov = curve_fit(func, x, ymean,[At, 2, 1], maxfev=10000)
		k = func(0,*popt)
		densmean[slice]=k
		if k > 800:
			start = slice - i
			i = i + 1

	end = start + i
	densmean2 = np.zeros(i + 1)
	i2 = 0
	for slice2 in range(0,200):
		ymean = input[slice2]
		At = guess(ymean)
		popt, pcov = curve_fit(func, x, ymean,[At, 2, 1], maxfev=10000)
		k = func(0,*popt)
		if k > 800:
			densmean2[i2] = k
			i2 = i2 + 1
	return np.median(densmean2)



# The function "checkslice" makes a fit of the density data, and the checks if the value of 
# the fitted function at R=0 is greater then 80% of the bulk value. If so, it returns "True". If not, it returns "False".
# "maxfev" is the allowed number of funtion calls when doing the fit. The default number was 800.
# It takes 3 arguments as input:
# s: (type: list) density data of a specific slice i, in our program it would be "input[i]"
# t: (type: float) the initial guess of the first parameter of the fit.
def checkslice(s, t, bulk):
	popt, pcov = curve_fit(func, x, s,[t, 2, 1], maxfev=10000)
	k = func(0,*popt)
	bulklimit = 0.9*bulk
	if k > bulklimit:
		return True
	else: 
		return False


# FUNCTION "from_slice_circle" makes a circle fit of some consecutive slices
# slicemin: number of the first slice
# slicemax: number of the last slice
# input: array containing all the data of the density map of all the slices 
# zcoord: list with the z-coordinates of all the slices
# returns 7 variables: 
# xc, yc, R_2: the coordinates of the circle center and its radius (xc corresponds to the z-position and yc to the radial distance to the central axis of the drop). 
# But yc won't be used for later calculations and has to be set to yc=0).
# residu, residu2: error and quadratic error of the circle fit 
# z, R: coordinates of the resulting points of the fitted radial density profile of each slice between slicemin and slicemax
# We use part of an imported source code to fit points to a circle. Source code can be found here--> "http://wiki.scipy.org/Cookbook/Least_Squares_Circle"

#the functions calc_R and f_3 will be used in the function from_slice_circle
def calc_R(xc, yc):
	""" calculate the distance of each 2D points from the center c=(xc, yc) """
	return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_3(beta, x):
	""" implicit definition of the circle """
	return (x[0]-beta[0])**2 + (x[1])**2 -beta[2]**2
	
def from_slice_circle(slicemax, slicemin, input, zcoord):
# Now we determine R and z. 
# R comes from the fit, and z is a vector with the shifted positions of the slices.
	dims = slicemax - slicemin + 1
	if dims<0:
		return
	else:
	
	 	x = input[0]
	 	i=0
	 	z = np.zeros(dims)
	 	R = np.zeros(dims)
	 	for a in range(slicemin, slicemax +1):
	 		y2 = input[a]
	 		i = a - slicemin
	 		z[i] = zcoord[a-1] 
	 		At = guess(y2)
	 		params, pcov2 = curve_fit(func, x, y2,[At, 2, 1], maxfev=10000)
	 		R[i] = params[1]
	 	
	 	# coordinates of the barycenter
	 	x_m = np.mean(z)
	 	y_m = np.mean(R)
	 	x = z
	 	y = R
	 
	 	# initial guess for parameters
	 	R_m = calc_R(x_m, 0).mean()
	 	beta0 = [ x_m, 0, R_m]
	 
	 	# for implicit function :
	 	# data.x contains both coordinates of the points
	 	# data.y is the dimensionality of the response
	 	lsc_data= odr.Data(np.row_stack([x, y]), y=1)
	 	lsc_model  = odr.Model(f_3, implicit=True)
	 	lsc_odr	   = odr.ODR(lsc_data, lsc_model, beta0)
	 	lsc_out	   = lsc_odr.run()
	 
	 	xc_3, yc_3, R_3 = lsc_out.beta
	 	Ri_3 = calc_R(xc_3, yc_3)
	 	residu_3 = sum((Ri_3 - R_3)**2)
	 	residu2_3 = sum((Ri_3**2-R_3**2)**2)

	 	return xc_3, yc_3, R_3, residu_3, residu2_3, z, R



# The function "distance" calculates the distance between the fitted circle and a point. It takes 5 argumets:
# x, y are the coordinates of the point; xc, yc are the coordinates of the circle center; R is its radius.
def distance(x, y, xc, yc,R):
	a = np.asarray([x,y])
	b = np.asarray([xc,yc])
	d = np.linalg.norm(a - b)- R
	d = math.fabs(d)
	return d 

# The function "plot_circ" makes a plot of the fitted circle together with the points that have been fitted. 
# The arguments: z is the shifted position of each slice, R the calculated radius corresponding to each slice,
# xc and yc are the coordinates of the center of the fitted circle, and R_2 is the radius of the fitted circle.
# It returns the figure object "fig".
# Other options regarding the aspect of the plots should be added/changed here inside de function.
def plot_circ(z, R, xc, yc, R_2):
	fig = pl.figure()
	theta_fit = np.linspace(-np.pi, np.pi, 180)
	x_fit2 = xc + R_2*np.cos(theta_fit)
	y_fit2 = yc + R_2*np.sin(theta_fit)
	pl.plot(y_fit2, x_fit2, 'b-', R, z, 'k.')
	fig.suptitle(filename)
	pl.xlim(0,2*R_2+ 1.0)
	pl.ylim(0,2*R_2+ 1.0)
	pl.xlabel('R [nm]')
	pl.ylabel('z [nm]')
	pl.grid(True)
	return fig
	
# The function "angle" calculates the contact angle.
def angle(R, m):
	t = m / np.sqrt(R**2 - m**2)
	t2 = math.degrees(np.arctan(t)) + 90
	return t2
	
				
# The function "radius" calculates the base radius of the droplet.
def radius(R, m):
	r = np.sqrt(R**2 - m**2)
	return r					

