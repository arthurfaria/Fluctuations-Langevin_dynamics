
# coding: utf-8

# In[1]:


import numpy as np
import time as tm 
import pylab as plt

from scipy import integrate
from scipy import sparse


# In[2]:


################## Parameters value

dx = 1.0                     # spatial separation
x = np.arange(-50,50, dx)    # spatial grid points

m     = 1.0                  # mass
omega = 0.5                  # HO frequency
gamma = 10.0                 # dissip                        
Dif   = 2.0 * gamma/gamma    # D= nu/gamma^2 = 2KbT/gamma diffusion const
sigma = 2.0 * gamma          # range/gamma^2
kappa = gamma * Dif/(omega**2)

varsigma = sigma/gamma

beta0 = 1/2
beta1 = 1/8 * varsigma**2
beta2 = 1/48 * varsigma**4
beta3 = 1/384 * varsigma**6
beta4 = 1/3840 * varsigma**8
beta5 = 1/46080 * varsigma**10
beta6 = 1/645120 * varsigma**12

#print(Dif*beta0)
#print(Dif*beta1)
#print(Dif*beta2)
#print(Dif*beta3)
#print(Dif*beta4)
#print(Dif*beta5)
#print(Dif*beta6)

print("")

################## Initial rho0   

def ddf(x,sig):
    x0 = x #- 10
    val = np.zeros_like(x0)
    val[(-(1/(2*sig))<=x0) & (x0<=(1/(2*sig)))] = 1
    return val

delta0 = ddf(x,1)            # initial delta
print(np.trapz(delta0, x))


################## overdamped standard BM rho_t   

def BM_rho_t(t, rho):
    B0 = omega**2/gamma  *  np.convolve(x*rho, [1,-1], 'same') / dx  
    B2 = 1/2 *np.convolve(rho, [1,-2,1], 'same') / (dx**2)
    #print(t)
    return  B0 + Dif* B2 


################## overdamped GenBM rho_t   

def GenBM_rho_t(t, rho):
    D0 = omega**2/gamma  *  np.convolve(x*rho, [1,-1], 'same') / dx  
    D2 = 1/2 * np.convolve(rho, [1,-2,1], 'same') / (dx**2)
    D4 = 1/8 * varsigma**2 * np.convolve(rho, [1, -4, 6, -4, 1], 'same') / (dx**4)
    D6 = 1/48 * varsigma**4 * np.convolve(rho, [1, -6, 15, -20, 15,  -6, 1], 'same') / (dx**6)
    D8 = 1/384 * varsigma**6 * np.convolve(rho, [1, -8, 28, -56, 70, -56, 28, -8, 1], 'same') / (dx**8)
    D10 = 1/3840 * varsigma**8 * np.convolve(rho, [1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1], 'same') / (dx**10)
    return  D0 + Dif*(D2 + D4 + D6 + D8 + D10)

################## Harmonic Force and Potential    

V = 0.5 * omega**2 * (x)**2     # harmonic potential
#F = - omega**2 * (x)           # harmonic force

################## Plot of rho0 and V

#plt.plot(x, V, "k--", label=r"$V(x)$")           # plot potential
plt.plot(x, delta0, label=r"$\rho(x,0)$")         # plot rho0

#plt.xlim(-10, 10)
plt.legend(loc=1, fontsize=12, fancybox=False)
plt.figure()


# In[3]:


################## Time initial conditions

dt = 0.1                        # time interval for snapshots
t0 = 0.0                        # initial time
tf = 100.                       # final time
t_eval = np.arange(t0, tf, dt)  # time range

for i in t_eval:
    t_dur = int(i/dt)
    
#print(t_dur)

start = tm.time()

#########################################################

################## Solve the Initial Value Problem

# overdamped standard BM
BM_sol = integrate.solve_ivp(BM_rho_t, t_span=[t0, tf], y0=delta0, t_eval=t_eval,  method="RK23")

# overdamped generalized BM
GenBM_sol = integrate.solve_ivp(GenBM_rho_t, t_span=[t0, tf], y0=delta0, t_eval=t_eval,  method="RK23")

print(GenBM_sol.message)

#########################################################

end = tm.time()
print("")
print(end-start)


# In[4]:


################## Plot of BM num and BM an

fig1 = plt.figure(figsize=(8, 4))
kappa = gamma * Dif/omega**2

x0 = - x/(kappa)
Z = 1./(np.sqrt(2*np.pi*kappa * np.abs(- beta0)))
brow_an = Z * (np.exp(-(kappa)*x0**2))

#########

plt.plot(x, brow_an , color='red', label=r"$ an $")
plt.plot(x, BM_sol.y[:,t_dur], color= 'darkblue', label=r"$ num $")

plt.xlabel(r' $x$', fontsize = 12)
plt.ylabel(r' $\rho_{BM}(x)$', fontsize = 12)

plt.legend(loc=1, fontsize=12, fancybox=False)
#fig1.savefig('/home/fariaart/Downloads/compar.png', transparent=False)

plt.figure()


# In[5]:


################## Plot of GenBM and Brown rhos

fig2 = plt.figure(figsize=(8, 4))

######### Snapshots

#for i, t in enumerate(GenBM_sol.t):
    #if i != 0:
        #print(np.trapz(sol.y[:,i],x))
        #if t//2 == 1:
            #plt.plot(x, GenBM_sol.y[:,i])
            
        
plt.plot(x, GenBM_sol.y[:,t_dur], color='red', label=r"$\rho_{GenBM}(x)$")
plt.plot(x, BM_sol.y[:,t_dur], color= 'darkblue', label=r"$\rho_{BM}(x)$")

#plt.ylim(0,0.08)
#plt.xlim(-20, 20)
plt.xlabel(r' $x$', fontsize = 12)

plt.legend(loc=1, fontsize=12, fancybox=False)
#fig2.savefig('/home/fariaart/Dropbox/Pesquisa/Doutorado/Lutz/Steady/Sigma/rho_ss_%s.png' %varsigma, transparent=False)


# In[6]:


################## Moments

BM_m = np.trapz(x*BM_sol.y[:,t_dur],x)
BM_v = np.trapz(x**2*BM_sol.y[:,t_dur],x)
BM_k = (np.trapz(x**4*BM_sol.y[:,t_dur],x) - 3*BM_v**2)/BM_v**2

GenBM_m = np.trapz(x*GenBM_sol.y[:,t_dur],x)
GenBM_v = np.trapz(x**2*GenBM_sol.y[:,t_dur],x) 
GenBM_k = (np.trapz(x**4*GenBM_sol.y[:,t_dur],x) - 3*GenBM_v**2)/ GenBM_v**2


# Normalization

print(np.trapz(BM_sol.y[:,t_dur],x))
print(np.trapz(GenBM_sol.y[:,t_dur],x))
print("")

# Mean value

print(BM_m)
print(GenBM_m)
print("")

# Variance

print(BM_v)
print(GenBM_v)
print("")

# Kurtosis

print(BM_k)
print(GenBM_k)
print("")

