#find peak in graphene data

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import leastsq



    

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

lx = np.array([np.linspace(kx[0]-50,kx[0]+50,101)])
lx = np.append(lx,[np.linspace(kx[1]-50,kx[1]+50,101)],axis=0)
lx = np.append(lx,[np.linspace(kx[2]-50,kx[2]+50,101)],axis=0)          

k = np.array([0])
k = np.delete(k,0)

for j in lx:
    for i in j :
        k = np.append(k,gy[i])
        
ly = np.array([k[0:101]])
ly = np.append(ly, [k[101:202]], axis = 0)
ly = np.append(ly, [k[202:303]], axis = 0)

    
x_bg = np.array([0])
x_bg = np.delete(x_bg,0)
y_bg = np.array([0])
y_bg = np.delete(y_bg,0)

for j in lx : 
    ind_bg_low = (j > min(j)) & (j < j[19])
    ind_bg_high = (j > j[70]) & (j < max(j))

    x_bg = np.append(x_bg,np.concatenate((j[ind_bg_low],j[ind_bg_high])))
    

 
for j in x_bg :
    y_bg = np.append(y_bg,gy[j])
    
m, c = np.polyfit(x_bg, y_bg, 1)

# removing fitted background # 
background = m*lx + c
y_bg_corr = ly - background
#pylab.plot(x,y_bg_corr)

#########################################################################
############################# FITTING DATA ## ###########################

# initial values #
p = np.array([[20.0,kx[0],1e20]])
p = np.append(p,[[20.0,kx[1],1e20]],axis=0)
p = np.append(p,[[20.0,kx[2],1e20]],axis=0)  # [hwhm, peak center, intensity] #

def lorentzian(x,p):
    numerator =  (p[0]**2 )
    denominator = (x-(p[1]))**2 + p[0]**2
    y = p[2]*(numerator/denominator)
    return y

def residuals(p,y,x):
    err =  - lorentzian(x,p)
    return err
# optimization # 

ly1 = ly[0]
ly2 = ly[1]
ly3 = ly[2]

p1 = p[0]
p2 = p[1]
p3 = p[2]

lx1 = lx[0]
lx2 = lx[1]
lx3 = lx[2]

y_bg_corr1 = y_bg_corr[0]
y_bg_corr2 = y_bg_corr[1]
y_bg_corr3 = y_bg_corr[2]

lorentz1 = lorentzian(lx1,p1)
lorentz2 = lorentzian(lx2,p2) 
lorentz3 = lorentzian(lx3,p3) 

pbest1=leastsq(residuals,p1,args=(y_bg_corr1,lx1),full_output=1)
pbest2=leastsq(residuals,p2,args=(y_bg_corr2,lx2),full_output=1)
pbest3=leastsq(residuals,p3,args=(y_bg_corr3,lx3),full_output=1)

pbest = np.array([pbest1])
pbest = np.append(pbest,[pbest2],axis=0)
pbest = np.append(pbest,[pbest3],axis=0)

best_parameters1 = pbest[0][0]
best_parameters2 = pbest[1][0]
best_parameters3 = pbest[2][0]

# fit to data #
    
fit1 = -lorentzian(lx1,best_parameters1)
fit2 = lorentzian(lx2,best_parameters2)
fit3 = lorentzian(lx3,best_parameters3)     
fit = np.array([fit1])
fit = np.append(fit,[fit2])
fit = np.append(fit,[fit3])

lx23 = np.arange(lx2[-1],lx3[0],1)

nx = np.append(lx1,lx2)

nx = np.append(nx,lx3)

y_bg_c23 = np.linspace(0,0,num=(lx3[0]-lx2[-1]))
y_bg_c = np.append(y_bg_corr1,y_bg_corr2)
y_bg_c = np.append(y_bg_c,y_bg_c23)
y_bg_c = np.append(y_bg_c,y_bg_corr3)


plt.plot(lx,y_bg_corr,'wo')
plt.plot(nx,fit,'r-',lw=2)
plt.plot(gx,gy)
plt.plot(peakind,gy[peakind])
plt.plot(kx,ky)
plt.show()