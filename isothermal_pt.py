# This is an Isothermal PT model
import numpy as np
import matplotlib.pyplot as plt

V_c0 = 0.0001
V_cs = 0.1*V_c0
V_e0 = 0.00005
V_es = 0.1*V_e0
V_b0 = 0.00005
V_bs = 0.1*V_b0
phi = 60
pi = np.pi
f = 50
T_h = 300
T_c = 77
P_0 = 2*10**6

V_reg = 0.25*pi*0.045**2*0.065*0.7
V_hot = 0.1*V_reg
V_pul = 0.025*pi*0.024**2*0.1
V_cold = 0.1*V_pul
V_pon = 1.25
V_g0 = V_pul*V_pon
T_reg = (T_h - T_c)/np.log(T_h/T_c)
R = 2078.5
V_T = V_hot/T_h + V_reg/T_reg + V_cold/T_c
nt = 360
dt = 1/f/nt
nx = 26

n = np.linspace(1,nt,nt)
V_c = np.zeros(nt)
V_e = np.zeros(nt)
V_b = np.zeros(nt)
V_c[:] = V_cs + 0.5*V_c0*(1-np.cos(n[:]*pi/180-phi*pi/180))
V_e[:] = V_es + 0.5*V_e0*(1-np.cos(n[:]*pi/180))
V_b[:] = V_bs + 0.5*V_b0*(1+np.cos(n[:]*pi/180))

#--------------------
P = np.zeros((nt,nx))
P[:,0] = 1/((V_c+V_b)/T_h+V_T+(V_pul+V_e)/T_c)
P[:,1] = P[:,0]*P_0/np.mean(P[:,0])
for i in range(2,11,2):
    P[:,i] = 1/((V_c+V_b)/T_h+V_T+(V_pul+V_e-V_g0*(P[:,i-1]/P_0)**-0.6)/T_c)
    P[:,i+1] = P[:,i]*P_0/np.mean(P[:,i])

for i in range(12,25,2):
    P[:,i] = 1/((V_c+V_b)/T_h+V_T+(V_pul+V_e-V_g0*(0.5*(P[:,i-1]+P[:,i-5]))**-0.6)/T_c)
    P[:,i+1] = P[:,i]*P_0/np.mean(P[:,i])

V_ge = np.zeros(nt)
V_ge = V_pul + V_e - V_g0*(P[:,25]/P_0)**-0.6
M_cold = P[:,25]*(V_c+V_ge)/(R*T_c)
M_reg = M_cold + P[:,25]*V_reg/(R*T_reg)
M_frc = np.zeros(nt)
M_frh = np.zeros(nt)
deltaV_c = np.zeros(nt)
deltaV_e = np.zeros(nt)
for i in range(nt):
    M_frc[i] = (M_cold[i]-M_cold[i-1])/dt
    M_frh[i] = (M_reg[i]-M_reg[i-1])/dt
    deltaV_c[i] = V_c[i]-V_c[i-1]
    deltaV_e[i] = V_e[i]-V_c[i-1]

PV_c_sum = sum(P[:,25]*deltaV_c)
PV_e_sum = sum(P[:,25]*deltaV_e)

plt.figure(1)
l_Vc = plt.plot(n,V_c,label = 'V_c')
l_Ve = plt.plot(n,V_e,label = 'V_e')
l_Vge = plt.plot(n,V_ge,label = 'V_ge')
plt.legend()
plt.show()

plt.figure(2)
fig, ax1 = plt.subplots()
ax1.plot(n,M_frc,label = 'M_frc')
ax1.plot(n,M_frh,label = 'M_frh')
plt.legend()
ax2 = ax1.twinx()
ax2.plot(n,P[:,25],'r-.',label = 'P25')
plt.legend()
fig.tight_layout()
plt.show()

plt.figure(3)
plt.plot(V_c,P[:,25],label = 'V_c')
plt.plot(V_e,P[:,25],label = 'V_e')
plt.plot(V_ge,P[:,25],label = 'V_ge')
plt.legend()
plt.show()

deltaP = P[:,25] - P[:,23]
