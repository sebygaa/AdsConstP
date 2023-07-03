
#%%
# Importing Libraries
import numpy as np
from scipy.integrate import odeint

# %%
# Basic conditions
L = 1
D = 0.2
N = 101
A_cros = D**2/4*np.pi
epsi = 0.4      # m^3/m^3
rho_s = 1000    # kg/m^3
dp = 0.02       # m (particle diameter)
#mu = [1.81, 1.81] # Pa sec (visocisty of gas)
mu_av = 1.81E-5
n_comp = 2

# %%
# Pressure boundary conditions
P_out = 1       # bar
#### HOW TO DETERMINE THIS? ####
T_feed = 300    # K
P_feed = 2      # bar
#################################
u_feed = 0.1            # m/s
y_feed = [0.2, 0.8]

# %%
## Ergun equation
Q = u_feed*A_cros*epsi  # volumetric flowrate

Cv = 5E-3       # m^3/s/bar
P_L = P_out + Q/Cv
Rgas = 8.314    # J/mol/K
P_av = (P_out+P_feed)/2
DPDz_term1 = -1.75*P_av*1E5/Rgas/T_feed*(1-epsi)/epsi*u_feed**2
DPDz_term2 = -150*mu_av/dp**2*(1-epsi)**2/epsi**2*u_feed
DPDz = DPDz_term1 + DPDz_term2
print(P_av)
print(DPDz)
P_0 = P_L + DPDz*1E-5
P_in = P_0 + Q/Cv
print('P_L = ', P_L)
print('P_0 = ', P_0)
print('P_in = ', P_in)

# %%
## z domain and P for position
z = np.linspace(0,L, N)
P = DPDz*1E-5*z + P_0   # bar
C = P*1E5/Rgas/T_feed   # mol/m^3
print(P)
print (C)

# %%
# backward FDM
#h = L/(N-1)
h = []
for ii in range(N-1):
    h.append(z[ii+1]-z[ii])
h.append(h[-1])
#h_arr = h*np.ones(N)
h_arr = np.array(h)
    
d_00 = np.diag(1/h_arr)
d_lo = np.diag(-1/h_arr[:-1],-1)
d = d_00+ d_lo
d[0,0] = 0
print(d)

dd_00 = np.diag(-2/h_arr**2)
dd_hi = np.diag(1/h_arr[1:]**2, 1)
dd_lo = np.diag(1/h_arr[:-1]**2, -1)
dd = dd_00 + dd_hi + dd_lo
dd[0,0:2] = 0
dd[-1,-3:] = 0
dd[-1,-2] = 1/h_arr[-1]**2
dd[-1,-1] = -1/h_arr[-1]**2
print(dd)
# %% Isotherm function as T ad P
def isomix(P_in,T):
    pp = []
    for ii in range(len(P_in)):
        p_tmp = P_in[ii][:]
        #ind_tmp = P_in[ii] < 1E-4
        #p_tmp[ind_tmp] = 0
        pp.append(p_tmp)
    q1 = 1*pp[0]*0.05/(1 + 0.3*pp[0] + 0.05*pp[1])
    q2 = 3*pp[1]*0.3/(1 + 0.3*pp[0] + 0.05*pp[1])
    return [q1, q2]


# %%
# Info for mass transfer
k = [0.5, 0.5]
D_AB = [1E-6, 1E-6]
# %%
# ODE for mass balance
def massbal(y,t):
    y_comp = []
    P_part = []
    dy_comp = []
    ddy_comp = []
    for ii in range(n_comp-1):
        y_tmp = np.array(y[ii*N:(ii+1)*N])
        y_tmp[y_tmp<1E-6] = 0
        dy_tmp = d@y_tmp
        ddy_tmp = dd@y_tmp
        P_part_tmp = y_tmp*P
        #y_tmp[y_tmp<1E-6] = 0
        y_comp.append(y_tmp)
        P_part.append(P_part_tmp)
        dy_comp.append(dy_tmp)
        ddy_comp.append(ddy_tmp)

    y_rest = 1-np.sum(y_comp, 0)
    P_part.append(P*y_rest)
    #print('Shape of P_part[0]',P_part[0].shape)
    q = []
    for ii in range(n_comp-1,2*n_comp-1):
        q.append(y[ii*N:(ii+1)*N])
    # Solid uptake component
    dqdt = []
    qstar = isomix(P_part,300)
    for ii in range(n_comp):
        dqdt_tmp = k[ii]*(qstar[ii] - q[ii])
        dqdt.append(dqdt_tmp)
    # Solid uptake total
    Sig_dqdt = np.sum(dqdt, 0)
    #print(Sig_dqdt)
    dy_compdt = []
    for ii in range(n_comp-1):
        term1 = -u_feed*dy_comp[ii]
        term2 = D_AB[ii]*0 # How to calculate d(yC)dz
        term3 = -rho_s*(1-epsi)/epsi*dqdt[ii]
        term4 = 0
        term5 = +y_comp[ii]*rho_s*Sig_dqdt*(1-epsi)/epsi
        
        #term1[0] = 0
        term1[0] = u_feed*d[1,1]*(y_feed[ii]-y_comp[ii][0])
        #term5[0] = 0
        #term1[-1] = -u_feed*d[1,1]*C[-1]*(y_comp[ii][-2]-y_comp[ii][-1])
        #term3[-1] = 0
        #term5[-1] = 0

        #dydt_tmp = term1+(term2 + term3+ term4 + term5)/C
        dydt_tmp = term1 + term3/C + term5/C
        #dydt_tmp[0] = 0
        #dydt_tmp[N-1] = 0 
        #ind_both = ((y_comp[ii]<1E-6) + (dydt_tmp < 0))>1.5
        #dydt_tmp[ind_both] = 0 
        dy_compdt.append(dydt_tmp)
    
    dydt = []
    for yy in dy_compdt:
        dydt.append(yy)
        #print(yy.shape)
    for qq in dqdt:
        dydt.append(qq)
        #print(qq.shape)
    dydt_arr = np.concatenate(dydt)
    
    return dydt_arr

# %% 
# Initial value
y0 = np.zeros(N*3)
y0[:N] = 1.0
P0_1 = P_0*y0[0:N]
P0_2 = P_0*(1-y0[0:N])
print('Here')
q0_1,q0_2 = isomix([P0_1,P0_2], 300)

y0[N:2*N] = q0_1
y0[2*N:3*N] = q0_2

#y0[0]  = 0
massbal(y0, 1)
# %%
# Run
t_dom = np.linspace(0,100,10001)
y_res = odeint(massbal, y0,t_dom,)

# %%
print(y_res)
print(y_res.shape)
# %%

import matplotlib.pyplot as plt
ls = ['-','--','-.',':']
n_ls = -1
for ii in range(len(t_dom))[:5000:200]:
    n_ls = n_ls+1

    plt.plot(z, y_res[ii,0*N:1*N],)

# %%
n_ls = -1
for ii in range(len(t_dom))[:5000:200]:
    n_ls = n_ls+1

    plt.plot(z, y_res[ii,2*N:3*N],)

# %%
# HERE is an Important tip!
# Starting with high concentration makes easy to calculate
# Make the initial y1 high !
# %%