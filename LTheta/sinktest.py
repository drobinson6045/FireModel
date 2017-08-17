import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp
from matplotlib import pyplot as plt


def f(x):
    return -2*x**2 + 1

def fp(x):
    return -4*x

def integrand(x):
    return(np.sqrt(1+fp(x)**2))
global N
N = 63
xs = np.linspace(-1,1,N)
ys = f(xs) 
L = integrate.quad(integrand,-1,1)
print(L)
dx = 2.0/(N-1)
fLen = np.zeros(N)
total = 0
#Trapezoidal integration
for i in range(1,N):
    total += (integrand(xs[i])+integrand(xs[i-1]))/2*dx
    fLen[i] = total
lInterp = interp.pchip(fLen,xs)
xsN = np.zeros(N)
lengths = np.linspace(0,fLen[-1],N)
xsN = lInterp(lengths)
xs = xsN
ys = f(xsN)
[xp,yp] = np.meshgrid(np.linspace(-2,2,130),np.linspace(-2,2,130))
sinkx = np.zeros((len(xp),len(xp)))
sinky = np.zeros((len(xp),len(xp)))
print(xp.shape)
for i in range(len(xp)):
    for k in range(len(xp)):
        for p in range(N):
            r = np.sqrt((xp[i,k]-xs[p])**2 + (yp[i,k]-ys[p])**2)
            sinkx[i,k] -= (xp[i,k]-xs[p])/r**2
            sinky[i,k] -= (yp[i,k]-ys[p])/r**2
print(sinky[:,15])
sinky +=80.0
plt.streamplot(xp,yp,sinkx,sinky)
plt.plot(xs,ys)
#plt.scatter(xs,ys)
plt.axis('equal')
plt.show()








