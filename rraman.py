#! /usr/bin/env python

import numpy as np
import pylab 
from scipy.optimize import leastsq

def lorentzian(x,p):
    numerator =  (p[0]**2 )
    denominator = ( x - (p[1]) )**2 + p[0]**2
    y = p[2]*(numerator/denominator)
    return y

def gaussian(x,p):
    c = p[0] / 2 / np.sqrt(2*np.log(2))
    numerator = (x-p[1])**2
    denominator = 2*c**2
    y = p[2]*np.exp(-numerator/denominator)
    return y

def residuals(p,y,x):
    err = y - lorentzian(x,p)
    return err

def local_fit(x, y, section):
    x=x[section[0]:section[1]]
    y=y[section[0]:section[1]]
    y_bg=y.min()
    p = [(x.max()-x.min())/2, (x.max()-x.min())/2+x.min() , x.max()-x.min()]
    # [fwhm, peak center, intensity] #
    pbest = leastsq(residuals, p, args=(y-y_bg,x), full_output=1)
    best_parameters = pbest[0]
    best_parameters[0] *= 2 
    fit = lorentzian(x,best_parameters) + y_bg
    return best_parameters,  x, fit

    

def main(f_object):
    ys = np.loadtxt(f_object)
    x_orig = np.linspace(1240.0691, 2919.5986, 1024)

    sections = [[0, 150], [150, 300], [800, 1000]]
    fit_results = []
    for line in ys:
        fit_result_for_a_line = [line[0]]
        for section in sections:
            fit_result = local_fit(x_orig, line[1:], section)
            fit_result_for_a_line = np.append(fit_result_for_a_line, np.abs(fit_result[0]), axis = 0)
        fit_result_for_a_line = np.append(fit_result_for_a_line, fit_result_for_a_line[9] / fit_result_for_a_line[6])
        fit_results = np.append(fit_results, fit_result_for_a_line, axis=1)
        print(line[0])
    fit_results=np.reshape(fit_results,(ys[:,0].size,11))

    Header_text = 'Time\t' + 'D_FWHM\tD_Center\tD_Amplitude\t' + \
    'G_FWHM\tG_Center\tG_Amplitude\t' + '2D_FWHM\t2D_Center\t2D_Amplitude\t2D/G Ratio'
    np.savetxt(f_object+'.fitted.txt', fit_results, delimiter='\t', newline='\r\n', header=Header_text)
    
    line = ys[300,:]
    pylab.plot(x_orig, line[1:], 'wo')
    for section in sections:
        fit_result = local_fit(x_orig, line[1:], section)
        print(fit_result[0])
        pylab.plot(fit_result[1], fit_result[2], 'r-', lw=2)
    pylab.show(block=False)

if __name__ == '__main__':
    f_object = 'C:\\Users\\sec\\Desktop\\manual\\graphene_3-3-C_200slit_0.2sec_5.1desity'
    main(f_object)