# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp





#basic reproduction number
R_0 = 2.6

#initial dubbeling time
d_0 = 3.5

#initial growth rate
r_0 = np.log(2)/d_0

#mean generation time
g = (R_0-1)*d_0/np.log(2) 

#Current doubeling time
d = 9

#Effective growth rate
r =np.log(2)/d

#Effective contact rate
l = r + 1/g

#infection fatility risk
f = 0.02

#Social distancing measures
rho = (r_0-r)*R_0/g

#Effective reproduction number
R = 1+r*g 

#Population of Stockholm
N = int(2e6)

# recovered fraction 3 weeks ago
I = 0.03 

# factor of unreported deaths
D_u = 1

#Initial data 
i_0 = 1000*r/f 
r_0 = 1000*(1-r)/f
s_0 = N-i_0 -r_0

t0 = -21
tf = 200

t_span = (t0,tf)

def SIR(t,x):
    s, i, r = x[0], x[1]*( x[1] > 1e-5 ) , x[2]
    s_prime = - s * l * i/N
    i_prime = ( s * l/N - 1/g ) * i
    r_prime = 1/g*i
    return np.array([s_prime, i_prime, r_prime])


x_0 = np.array([s_0, i_0, r_0])

Sol = solve_ivp(SIR,t_span, x_0,method = 'BDF', dense_output =True)


Total = N-int(Sol.sol(tf)[0])
D = int(Total*f)

t = np.arange(t0,tf,1e-3)
t_d = np.arange(t0+21,tf,1e-3)
#plt.plot(t,Sol.sol(t)[0],label = ('susceptible') )

plt.plot(t,Sol.sol(t)[1],label = 'f={}, D_tot = {}'.format(f,D)  )
#plt.plot(t_d,Sol.sol(t_d-21)[1]*f, label = 'daily deaths' )
#plt.plot(t,Sol.sol(t)[2],label = ('removed') )
plt.xlabel('days from today')
plt.ylabel('infected')
plt.title('Stockholm, R_eff={}, d_eff ={} '.format(
    int(10*R)/10, d) )
plt.grid()
plt.legend(loc='upper right')
plt.show()
