import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz as w
p = np.pi
print('p=',p)

# Brendel-Bormann (BB) model parameters
wp = 14.98 #eV
f0 = 0.526 
L0 = 0.047 #eV

f1 = 0.213
L1 = 0.312 #eV
w1 = 0.163 #eV
o1 = 0.013 #eV

f2 = 0.060
L2 = 0.315 #eV
w2 = 1.561 #eV
o2 = 0.042 #eV

f3 = 0.182
L3 = 1.587 #eV
w3 = 1.827 #eV
o3 = 0.256 #eV

f4 = 0.014
L4 = 2.145 #eV
w4 = 4.495 #eV
o4 = 1.735 #eV

Omp = f0**.5 * wp  #eV

def BB(wl):  #wl: eV
    e = 1-Omp**2/(wl*(wl+1j*L0))

    a = (wl**2+1j*wl*L1)**.5
    za = (a-w1)/(2**.5*o1)
    zb = (a+w1)/(2**.5*o1)
    e += 1j*p**.5*f1*wp**2 / (2**1.5*a*o1) * (w(za)+w(zb)) #x1
    
    a = (wl**2+1j*wl*L2)**.5
    za = (a-w2)/(2**.5*o2)
    zb = (a+w2)/(2**.5*o2)
    e += 1j*p**.5*f2*wp**2 / (2**1.5*a*o2) * (w(za)+w(zb)) #x2
    
    a = (wl**2+1j*wl*L3)**.5
    za = (a-w3)/(2**.5*o3)
    zb = (a+w3)/(2**.5*o3)
    e += 1j*p**.5*f3*wp**2 / (2**1.5*a*o3) * (w(za)+w(zb)) #x3
    
    a = (wl**2+1j*wl*L4)**.5
    za = (a-w4)/(2**.5*o4)
    zb = (a+w4)/(2**.5*o4)
    e += 1j*p**.5*f4*wp**2 / (2**1.5*a*o4) * (w(za)+w(zb)) #x4
    
    return e
  
ev_min=0.005
ev_max=20
npoints=2000
eV = np.logspace(np.log10(ev_min), np.log10(ev_max), npoints)
um = 4.13566733e-1*2.99792458/eV
e = BB(eV)
n = (e**.5).real
k = (e**.5).imag


#============================   DATA OUTPUT   =================================
#file = open('out.txt', 'w')
#for i in range(npoints-1, -1, -1):
 #   file.write('\n        {:.4e} {:.4e} {:.4e}'.format(um[i],n[i],k[i]))
#file.close()
    
    
#===============================   PLOT   =====================================
plt.rc('font', family='Arial', size='14')

plt.figure(1)
plt.plot(eV, -e.real, label="-e1")
plt.plot(eV, e.imag, label="e2")
plt.xlabel('Photon energy (eV)')
plt.ylabel('e')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)


#plot n,k vs eV
plt.figure(2)
plt.plot(eV, n, label="n")
plt.plot(eV, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)


#plot n,k vs um
plt.figure(3)
plt.plot(um, n, label="n")
plt.plot(um, k, label="k")
plt.xlabel('Wavelength (um)')
plt.ylabel('n, k')

plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.show()