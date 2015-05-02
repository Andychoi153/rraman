# rraman and 3d dimenssion - plus
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
    
    for line in ys:
        
        for section in sections:
            fit_result = local_fit(x_orig, line[1:], section)
            print(fit_result[0])
            if section == [0,150] :
                fit_x = fit_result[1]
                fit_y = fit_result[2]
            elif section == [800,1000] :
                fit_x = np.append(fit_x, fit_result[1],axis=0)               
                fit_y = np.append(fit_y,np.linspace(245,245,500).tolist(),axis=0)
                fit_y = np.append(fit_y, fit_result[2],axis=0)
            else :
                fit_x = np.append(fit_x, fit_result[1],axis=0)
                fit_y = np.append(fit_y, fit_result[2],axis=0)
            
        if line[0] == ys[0][0]:
            totalfit_x = [fit_x]
            totalfit_y = [fit_y]
        else :
            totalfit_x = np.append(totalfit_x, [fit_x], axis=0)
            totalfit_y = np.append(totalfit_y, [fit_y], axis=0)
    return totalfit_x, totalfit_y
if __name__ == '__main__':
    import time
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    f_object = 'C:\\Users\\sec\\Desktop\\manual\\graphene_3-3-C_200slit_0.2sec_5.1desity'
    start = time.time()    
    x,z = main(f_object)
    end = time.time()
    i = range(1135)

            
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    bx = fig.gca(projection='3d')
    
    for i in range(1135):
        if i ==0 :
            X = np.array([np.linspace(1,1000,1000)])
        else :
            X = np.append(X,[np.linspace(1,1000,1000)],axis=0)
        
    
    for i in range(1135):
        if i ==0 :
            Y = np.array([np.linspace(1,1,1000)])
        else :
            Y = np.append(Y,[np.linspace(i+1,i+1,1000)],axis=0)

    X = X[300:]
    Y = Y[300:]
    Z = z[300:]
    ax.plot_surface(X, Y, Z, rstride=50, cstride=50, alpha=0.2)


    cset = ax.contourf(X, Y, Z, zdir='z', offset=220, cmap=cm.coolwarm)
    cset = ax.contourf(X, Y, Z, zdir='y', offset=1200, cmap=cm.coolwarm)

    ax.set_xlabel('X')
    ax.set_xlim(0, 1100)
    ax.set_ylabel('Y')
    ax.set_ylim(0,1200 )
    ax.set_zlabel('Z')
    ax.set_zlim(220, 300)

    plt.show()
    print (start-end)