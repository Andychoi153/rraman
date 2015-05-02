
#find peak in graphene data

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import leastsq


def lorentzian(lx,p):
    for j in p :
        numerator =  (j[0]**2 )
        denominator = ( lx - (j[1]) )**2 + j[0]**2
        ly = j[2]*(numerator/denominator)
    return ly

def residuals(p,ly,lx):
    err = ly - lorentzian(lx,p)
    return err
    

f= 'C:\\Users\\sec\\Desktop\\scan\\graphene_3-3-C_200slit_0.2sec_5.1desity'

x = np.loadtxt(f)
k = x[400]
gy = k[1:1024]
gx= np.linspace(0,1*gy.size, gy.size)
peakind = signal.find_peaks_cwt(gy,np.arange(1,80))

a = gy[peakind]
s = np.argsort(a)[-4:]
vals = a[s]
mx = np.array([0])

mx = np.delete(mx,0)
for i in s:
    mx = np.append(mx,peakind[i])

kx = np.array([0])
kx = np.delete(kx,0)
    
for j in mx :
    kx = np.append(kx,[(j-10)+np.argmax(gy[j-10:j+10])])
    
kx = kx[np.argsort(kx)]
kx = np.delete(kx,2)
ky = gy[kx]

lx = np.array([0])
lx = np.delete(lx,0)

for j in kx :
    lx = np.append(kx,[np.linspace(j-50,j+50,100)]) #lorentzian x


ly = np.array([0])
ly = np.delete(ly,0)

for j in lx:
    ly = np.append(ly,[gy[j]])

x_bg = np.array([0])
x_bg = np.delete(x_bg,0)
y_bg = np.array([0])
y_bg = np.delete(y_bg,0)

for j in lx : 
    ind_bg_low = (j > min(j)) & (j < j[19])
    ind_bg_high = (j > j[70]) & (j < max(j))

    x_bg = np.append(x_bg,np.concatenate((j[ind_bg_low],j[ind_bg_high])))
    
    
for j in x_bg :
    y_bg = np.append(y_bg,ly[j])
    
#pylab.plot(x_bg,y_bg)

# fitting the background to a line # 
m, c = np.polyfit(x_bg, y_bg, 1)

# removing fitted background # 
background = m*lx + c
y_bg_corr = ly - background
#pylab.plot(x,y_bg_corr)

#########################################################################
############################# FITTING DATA ## ###########################

# initial values #
p = np.array([0])
p = np.delete(p,0)
for j in kx :
    p = np.append(p,[[20.0,j,12e3]])  # [hwhm, peak center, intensity] #

# optimization # 
pbest = leastsq(residuals,p,args=(y_bg_corr,lx),full_output=1)

best_parameters = np.array([0])
best_parameters = np.delete(p,0)

for j in pbest :
    best_parameters = np.append(best_parameters, j[0])
# fit to data #
    
fit = lorentzian(lx,best_parameters)
    

plt.plot(gx,gy)
plt.plot(peakind,gy[peakind])
plt.plot(kx,ky)
plt.show()

