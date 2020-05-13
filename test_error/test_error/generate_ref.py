#from sympy.functions.special.spherical_harmonics import Znm
from sympy import assoc_legendre
from sympy import N, re,sin,cos,pi,factorial, sqrt, Pow
from sympy import symbols
from numpy import random
s_theta, s_phi = symbols('theta phi')
n_points = 1000
theta = random.rand(n_points) * pi
phi = random.rand(n_points) * 2 * pi
def real_sph(l,m,theta,phi):
    if m < 0:
        return Pow(-1,m)*sqrt(2)*sqrt((2*l+1)/(4*pi)*factorial(l+m)/factorial(l-m)) \
                *assoc_legendre(l,-m,cos(theta))*sin(-m*phi)
    if m == 0:
        return sqrt((2*l+1) / (4*pi)) * assoc_legendre(l,m,cos(theta))
    if m > 0:
        return Pow(-1,m)*sqrt(2)*sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m)) \
                *assoc_legendre(l,m,cos(theta))*cos(m*phi)

for i in range(n_points):
    print(N(theta[i],50),N(phi[i],50),end = ' ')
    print(N(sin(theta[i]),50),N(cos(theta[i]),50),N(sin(phi[i]),50),N(cos(phi[i]),50), end = " ")
    for l in range(21):
        for m in range(-l,l+1):
            #print(re(N(Znm(l,m,theta[i],phi[i]),50)),end = ' ')
            print(N(real_sph(l,m,s_theta,s_phi).subs([(s_theta,N(theta[i],50)),(s_phi,N(phi[i],50))]),50),end = ' ')
    print()
