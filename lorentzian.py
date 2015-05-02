#########################################################################
####################### IMPORTING REQUIRED MODULES ######################

import numpy
import pylab
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #

#########################################################################
############################# LOADING DATA ##############################

a = numpy.loadtxt('silicon_Raman.txt')
x = a[:,0]
y = a[:,1]

#########################################################################
########################### DEFINING FUNCTIONS ##########################

def lorentzian(x,p):
    numerator =  (p[0]**2 )
    denominator = ( x - (p[1]) )**2 + p[0]**2
    y = p[2]*(numerator/denominator)
    return y

def residuals(p,y,x):
    err = y - lorentzian(x,p)
    return err

#########################################################################
######################## BACKGROUND SUBTRACTION #########################

# defining the 'background' part of the spectrum #
ind_bg_low = (x > min(x)) & (x < 450.0)
ind_bg_high = (x > 590.0) & (x < max(x))

x_bg = numpy.concatenate((x[ind_bg_low],x[ind_bg_high]))
y_bg = numpy.concatenate((y[ind_bg_low],y[ind_bg_high]))
#pylab.plot(x_bg,y_bg)

# fitting the background to a line # 
m, c = numpy.polyfit(x_bg, y_bg, 1)

# removing fitted background # 
background = m*x + c
y_bg_corr = y - background
#pylab.plot(x,y_bg_corr)

#########################################################################
############################# FITTING DATA ## ###########################

# initial values #
p = [5.0,520.0,12e3]  # [hwhm, peak center, intensity] #

# optimization # 
pbest = leastsq(residuals,p,args=(y_bg_corr,x),full_output=1)
best_parameters = pbest[0]

# fit to data #
fit = lorentzian(x,best_parameters)

#########################################################################
############################## PLOTTING #################################

pylab.plot(x,y_bg_corr,'wo')
pylab.plot(x,fit,'r-',lw=2)
pylab.xlabel(r'$\omega$ (cm$^{-1}$)', fontsize=18)
pylab.ylabel('Intensity (a.u.)', fontsize=18)

pylab.show()