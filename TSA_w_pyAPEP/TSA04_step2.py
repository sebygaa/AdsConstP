# Determine Proper Cv values
# using the simulations with no mass transfer

from pyapep.simsep_rev import column
import numpy as np
import matplotlib.pyplot as plt
def Arr(T, dH, T_ref):
    thet = np.exp(dH/8.3145*(1/T - 1/T_ref))
    return thet
qm_N2 = 0.01 # mol/kg
b_N2 = 0.05 # bar^-1
dH_N2 = 750 # J/mol
#def iso_N2(P,T):
#    P_norm = P*Arr(T, dH_N2, 300)
#    numer = qm_N2 * b_N2*P_norm
#    denom = 1 + b_N2*P_norm
#    qN2 = numer/denom
#    return qN2

qm_O2 = 0.02 # mol/kg
b_O2 = 0.05 # bar^-1
dH_O2 = 850 # J/mol

qm_H2O = 1.5 # mol/kg
b_H2O = 0.2 # bar^-1
dH_H2O = 1600 # J/mol

qm_CO2 = 5.5 # mol/kg
b_CO2 = 2.7 # bar^-1
dH_CO2 = 1100 # J/mol

def iso_mix(P,T,):
    P1,P2,P3,P4 = P
    P1_norm = P1*Arr(T, dH_N2, 300)
    P2_norm = P2*Arr(T, dH_O2, 300)
    P3_norm = P3*Arr(T, dH_H2O, 300)
    P4_norm = P4*Arr(T, dH_CO2, 300)
    numer1 = qm_N2*b_N2*P1_norm
    numer2 = qm_O2*b_O2*P2_norm
    numer3 = qm_H2O*b_H2O*P3_norm
    numer4 = qm_CO2*b_CO2*P4_norm
    denom = 1 + b_N2*P1_norm + b_O2*P2_norm + b_H2O*P3_norm + b_CO2*P4_norm
    qN2 = numer1/denom
    qO2 = numer2/denom
    qH2O = numer3/denom
    qCO2 = numer4/denom
    return qN2, qO2, qH2O, qCO2

L = 4.0 # m
A_cr = (8**2)*1/360*np.pi # m^2

N = 41
c1 = column(L, A_cr, 4, N, )
c1.adsorbent_info(iso_mix, epsi = 0.4,
                  D_particle = 0.01, rho_s = 1000,)

Mw_arr = np.array([0.028, 0.032, 0.018, 0.044]) 
mu_arr = np.array([1.8E-5,]*4) # Pa.s
c1.gas_prop_info(Mw_arr,mu_arr)

k_arr = np.array([0.01,]*4)
c1.mass_trans_info(k_arr, a_specific_surf = 100,)

dH_arr = np.array([dH_N2, dH_O2, dH_H2O, dH_CO2])
Cp_s = 700 # J/(kg.K)
Cp_g = np.array([29.1, 29.4, 75.3, 37.1]) # J/(mol.K)
h_HTC = 0.001 # W/(m^2.K)
c1.thermal_info(dH_arr, Cp_s, Cp_g, h_HTC,)

y_inlet = np.array([0.7799, 0.20, 0.02, 0.0001])
Q_in = 112*1/270 # m^3/s modified
Cv_outlet = 55.0

T_feed = 298 # K

c1.boundaryC_info(P_outlet = 0.98, P_inlet = 1.6,
                   T_inlet = T_feed, y_inlet = y_inlet,
                   Q_inlet = Q_in, Cv_out = Cv_outlet)
P_init = 1.02*np.ones(N)
y_init = [0.7799*np.ones(N), 0.2001*np.ones(N),
          0.02*np.ones(N), 0.0000*np.ones(N),]
Tg_init = 298*np.ones(N)
Ts_init = 298*np.ones(N)
P_partial_init = [P_init*y_init[0], P_init*y_init[1],
                  P_init*y_init[2], P_init*y_init[3],]
q_init = iso_mix(P_partial_init, Tg_init)
c1.initialC_info(P_init, Tg_init, Ts_init,y_init, q_init)

# If the rotating period is 1000 seconds...
# Adsorption may take 750 seconds
y_res = c1.run_mamoen(750, 5, CPUtime_print = True)

#c1.Graph(10, 3, yaxis_label = 'CO$_{2}$',
#          file_name = 'TSA04_step1_C3.png', dpi = 300)
#c1.Graph(10, 3, yaxis_label = 'H$_{2}$O',
#         file_name = 'TSA04_step1_C4.png', dpi = 300)
t_evry = 50
c1.Graph(t_evry, 8, yaxis_label = 'T$_{g}$ (K)',
         file_name = 'TSA04_step1_Tg.png', dpi = 300)
c1.Graph(t_evry, 9, yaxis_label = 'T$_{s}$ (K)',
         file_name = 'TSA04_step1_Ts.png', dpi = 300)
#c1.Graph_P(10, file_name = 'TSA04_step1_P.png', dpi = 300)

c1.next_init(True)
# =-=-=-=-=-=-=-=-=-=-=-=-=
# Pressure stabilization
# =-=-=-=-=-=-=-=-=-=-=-=-=
y_inlet = np.array([0.01, 0.01, 0.98, 0.0])
Q_in = 0
Cv_outlet = 0
# =-=-=-=-=-=-=-=-=-=-=
Cv_inlet = 0
# =-=-=-=-=-=-=-=-=-=-=
T_feed = 373 # K


c1.boundaryC_info(P_outlet = 0.98, P_inlet = 1.6,
                   T_inlet = T_feed, y_inlet = y_inlet,
                   Q_inlet = Q_in, Cv_out = Cv_outlet, Cv_in = Cv_inlet,
                   foward_flow_direction= True)

c1.run_mamoen(2, 5, CPUtime_print = True)
#c1.Graph(0.5, 3, yaxis_label = 'CO$_{2}$',
#          file_name = 'TSA04_step2_C3.png', dpi = 300)
#c1.Graph(0.5, 3, yaxis_label = 'H$_{2}$O',
#         file_name = 'TSA04_step2_C4.png', dpi = 300)
#c1.Graph_P(1, file_name = 'TSA04_step2_P.png', dpi = 300)

# =-=-=-=-=-=-=-=-=-=-=-=-=
# Steam injection
# =-=-=-=-=-=-=-=-=-=-=-=-=
c1.next_init(True)

y_inlet = np.array([0.01, 0.01, 0.98, 0.0])
#Q_in_he = 112*1/270/2000 # m^3/s modified 1/2000 of flowrate ?
Cv_outlet_he = 7E-1
Cv_inlet_he = 7E-1
T_feed_he = 373 # K

k_arr_he = np.array([0.01,]*4)
c1.mass_trans_info(k_arr_he, a_specific_surf = 100,)

dH_arr = np.array([dH_N2, dH_O2, dH_H2O, dH_CO2])
Cp_s = 700 # J/(kg.K)
Cp_g = np.array([29.1, 29.4, 75.3, 37.1]) # J/(mol.K)
h_HTC = 1E1 # W/(m^2.K)
#h_HTC = 0.000 # W/(m^2.K)
c1.thermal_info(dH_arr, Cp_s, Cp_g, h_HTC,)

#c1.boundaryC_info(P_outlet = 0.98, P_inlet = 2.2,
#                   T_inlet = T_feed_he, y_inlet = y_inlet,
#                   Q_inlet = Q_in_he, Cv_out = Cv_outlet_he,
#                   foward_flow_direction= True)
c1.boundaryC_info(P_outlet = 0.98, P_inlet = 2.2,
                   T_inlet = T_feed_he, y_inlet = y_inlet,
                   Cv_out = Cv_outlet_he, Cv_in = Cv_inlet_he,
                   foward_flow_direction= False)
print(c1)

# If the rotating period is 1000 seconds...
# Regeneration step may take 250 seconds
t_evry = 25
c1.run_mamoen(t_evry*10-2, 5, CPUtime_print = True)

c1.Graph(t_evry, 3, yaxis_label = 'CO$_{2}$ conc. (mol/m$^{3}$)',
          file_name = 'TSA04_step3_C3.png', dpi = 300)
c1.Graph(t_evry, 2, yaxis_label = 'H$_{2}$O conc. (mol/m$^{3}$)',
         file_name = 'TSA04_step3_C2.png', dpi = 300)
c1.Graph(t_evry, 7, yaxis_label = 'CO$_{2}$ uptake (mol/kg)',
         file_name = 'TSA04_step3_q4.png', dpi = 300)
c1.Graph(t_evry, 8, yaxis_label = 'T$_{g}$ (K)',
         file_name = 'TSA04_step3_Tg.png', dpi = 300)
c1.Graph(t_evry, 9, yaxis_label = 'T$_{s}$ (K)',
         file_name = 'TSA04_step3_Ts.png', dpi = 300)
c1.Graph_P(t_evry, file_name = 'TSA04_step3_P.png', dpi = 300)
c1.Q_valve(True)

plt.show()


'''
Q_in_he = 112*1/270/100 # m^3/s modified 1/100 of flowrate
y_inlet_he = np.array([0.01, 0.01, 0.98, 0.0])
Cv_outlet_he = 350.0
T_feed_he = 398 # K
for ii in range(200):
    # =-=-=-=-=-=-=-=-=-=-=-=-=
    # Adsorption step
    # =-=-=-=-=-=-=-=-=-=-=-=-=
    
    # =-=-=-=-=-=-=-=-=-=-=-=-=
    # Heatting step
    # =-=-=-=-=-=-=-=-=-=-=-=-=
    
    c1.next_init(True)
    c1.boundaryC_info(P_outlet = 0.98, P_inlet = 2.3,
                       T_inlet = T_feed_he, y_inlet = y_inlet_he,
                       Q_inlet = Q_in_he, Cv_out = Cv_outlet_he,
                       foward_flow_direction= False)

    c1.run_mamoen(250-2, 5, CPUtime_print = True)
'''