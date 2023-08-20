# %%
# Importing the run_ma_linvel
from pyapep.simsep import *
import numpy as np
# %%
N = 101

D = 0.2
A_cros = D**2/4*np.pi
epsi = 0.4      # m^3/m^3
rho_s = 1000    # kg/m^3
dp = 0.02       # m (particle diameter)
mu = [1.81E-5, 1.81E-5] # Pa sec (visocisty of gas)
#mu_av = 1.81E-5
n_comp = 2

# %%
# Column Geometry
L = 1
D = 0.2
N = 41
A_cros = D**2/4*np.pi
epsi = 0.4      # m^3/m^3
# %%
# Define a column    
c1 = column(L, A_cros,n_comp, N, E_balance=False)

# %%
# Adsorbent info
def isomix(P_in,T):
    pp = []
    for ii in range(len(P_in)):
        p_tmp = P_in[ii][:]
        #ind_tmp = P_in[ii] < 1E-4
        #p_tmp[ind_tmp] = 0
        pp.append(p_tmp)
    q1 = 1*pp[0]*0.1/(1 + 0.1*pp[0] + 0.4*pp[1])
    q2 = 4*pp[1]*0.4/(1 + 0.1*pp[0] + 0.4*pp[1])
    return [q1, q2,]
c1.adsorbent_info(isomix, epsi, dp, rho_s,)

# %%
# Gas properties
Mw = [28, 32]
c1.gas_prop_info(Mw,mu)

# %%
# Mass transfer
k_MTC = [2.5, 2.5]
D_disp = [1E-7, 1E-7] 
a_surf = 1
c1.mass_trans_info(k_MTC, a_surf, D_disp)

# %%
# Boundary conditions
P_out = 1.0     # Ignore This Value
P_in = 1.0      # Important (during pressurization)
T_in = 300
y_in = [0.3, 0.7,]
Cv_in = 0        # m^3/sec/bar
Cv_out = 5E-2          # m^3/sec/bar
u_feed = 0.1            # m/s
Q_in = u_feed*A_cros*epsi  # volumetric flowrate
c1.boundaryC_info(P_out, P_in, T_in, y_in, Cv_in,Cv_out,Q_in, 
                    assigned_v_option = True, foward_flow_direction=False)

# %%
# Initial conditions
P_init = 5.5*np.ones([N,])

y_init = [0.3*np.ones([N,]),
            0.7*np.ones([N,]),]
P_part = [0.3*np.ones([N,])*P_init,
            0.7*np.ones([N,])*P_init,]
Tg_init = 300*np.ones([N,])
Ts_init = 300*np.ones([N,])
q1,q2 = isomix(P_part, Tg_init)
q_init = [q1,q2]

c1.initialC_info(P_init, Tg_init, Ts_init, y_init, q_init )
print(c1)

# %%
# Run
y_res, z_res, t_res = c1.run_ma_linvel(30,10)
# %%
c1.Graph_P(5)
# %%
# Graph of 1st component
c1.Graph(5, 3, loc = [1.13,0.8], 
         yaxis_label = 'Solid uptake of N$_2$ (mol/kg)')
# %%
c1.Graph(20, 1)
# %%
c1.Graph(20,2)
# %%
c1.Graph(10,5, )
# %%
plt.plot(c1._z, c1._y_fra[0][0,:N])
plt.plot(c1._z, c1._y_fra[0][1,:N])
plt.plot(c1._z, c1._y_fra[0][20,:N])
plt.plot(c1._z, c1._y_fra[1][0, :N])
plt.plot(c1._z, c1._y_fra[1][20, :N])
# %%
