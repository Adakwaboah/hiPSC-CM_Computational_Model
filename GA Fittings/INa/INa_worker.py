'''INa worker file - Luo-Rudy Model'''

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

dt = 0.02#time steps
tSpan = (0, 20)
tEval = np.arange(tSpan[0], tSpan[1], dt)

Vcm = np.arange(-80, 40, 5)
Vhold = -120  #mV


Na_o = 140 #mM

Na_i = 10  #mM
R = 8.314 #J K^-1 mol^-1
T = 310 #K
F = 96.4867 #C.mmol^-1

Ena = ((R*T)/F)*np.log(Na_o/Na_i)

#stimulus - aperiodic
#S1
S1_s = 10
S1_e = 35
S2_s = 110
S2_e = 135
#S3
S3_s = 210
S3_e = 235


Mthd = 'BDF'
# Mthd = 'LSODA'
# Mthd = 'DOP853'
# Mthd = 'RK45'
# Mthd = 'Radau'

frq = 20

#voltage clamping for stimulus
# def Vm(t, Vcm): 
#     if S1_s <= t < S1_e:
#         return Vcm
#     elif S2_s <= t < S2_e:
#         return Vcm
#     elif S3_s <= t < S3_e:
#         return Vcm
#     else:
#         return Vhold

def Vm(t, Vcm, freq):
    T = 1000/freq
    if abs(t)%T <= T/2:
        return Vcm
    else:
        return Vhold
    

#voltage stimulus generation - include frequency and duration - mainly for pacing hiPSC-CM have automaticity (so not needed anyway)
def Vstim(dur, Va, freq):
    Vstim = np.ones(len(dur))*Vhold
    T = 1000/freq
    i = 0
    for t in dur:
        if abs(t)%T <= T/2:
            Vstim[i] = Va
        i = i + 1
    return Vstim
    
m_o = 0
h_o = 1
j_o = 1
ini_con = (m_o, h_o, j_o)
             
def INa_dt(t, Y, Vc, INa_par):
    m = Y[0]
    h = Y[1]
    j = Y[2]

    am = INa_par[1]*(Vm(t, Vc, frq)+47.13)/(1-np.exp(-0.1*(Vm(t, Vc, frq)+47.13)))
    bm = INa_par[2]*np.exp(-Vm(t, Vc, frq)/11)
    
    if (Vm(t, Vc, frq) < -40):
        ah = INa_par[3]*np.exp((80+Vm(t, Vc, frq))/-6.8)
        bh = INa_par[4]*np.exp(0.079*(Vm(t, Vc, frq)))+INa_par[5]*np.exp(0.35*(Vm(t, Vc, frq)))
        aj = (-INa_par[6]*np.exp(INa_par[7]*Vm(t, Vc, frq))-INa_par[8]*np.exp(-INa_par[9]*Vm(t, Vc, frq)))*((Vm(t, Vc, frq)+37.78)/(1+np.exp(INa_par[10]*(Vm(t, Vc, frq)+79.23))))
        bj = (INa_par[11]*np.exp(-INa_par[12]*Vm(t, Vc, frq)))/(1+np.exp(-0.1378*(Vm(t, Vc, frq)+40.14)))
    else:
        ah = 0
        bh = 1/(INa_par[13]*(1+np.exp((Vm(t, Vc, frq)+INa_par[14])/-INa_par[15])))
        aj = 0
        bj = (INa_par[16]*np.exp(-INa_par[17]*Vm(t, Vc, frq)))/(1+np.exp(-INa_par[18]*(Vm(t, Vc, frq)+INa_par[19])))


    taum = 1/(am+bm)
    tauh = 1/(ah+bh)
    tauj = 1/(aj+bj)
    
    mss = am*taum
    hss = ah*tauh
    jss = aj*tauj
    
    dm_dt = ((mss-m)/taum)
    dh_dt = ((hss-h)/tauh)
    dj_dt = ((jss-j)/tauj)

    return dm_dt,dh_dt,dj_dt
# mhj_sol = solve_ivp(INa_dt, tSpan, ini_con, args = (Vc,) teval = tEval, method='Radau', dense_output=True)
# m = mhj_sol.y[0]
# h = mhj_sol.y[1]
# j = mhj_sol.y[2]
Cmem = 1
pts_limit = 10000
plot_INa = np.zeros((len(Vcm), pts_limit))
plot_m = np.zeros((len(Vcm), pts_limit))
plot_h = np.zeros((len(Vcm), pts_limit))
plot_j = np.zeros((len(Vcm), pts_limit))
time = np.zeros((len(Vcm), pts_limit))
INa_IVs = np.zeros(len(Vcm))
size_arr = np.arange(len(Vcm))

def mhj(parm):
    i = 0
    for Va in Vcm:
        mhj_sol = solve_ivp(INa_dt, tSpan, ini_con, args= (Va, parm), teval = tEval, method= Mthd)
#         Vstim = [Va if (S1_s <= k < S1_e) or (S2_s <= k < S2_e) or (S3_s <= k < S3_e) else Vhold for k in mhj_sol.t]
        Vmem = Vstim(mhj_sol.t, Va, frq) 
#         plt.plot(mhj_sol.t, Vmem)
        plot_m[i, :len(mhj_sol.y[0])] = mhj_sol.y[0]
        plot_h[i, :len(mhj_sol.y[1])] = mhj_sol.y[1]
        plot_j[i, :len(mhj_sol.y[2])] = mhj_sol.y[2]
        time[i, :len(mhj_sol.t)] = mhj_sol.t
        size_arr[i] = np.uint16(len(mhj_sol.t))
        i = i+1
    return time, size_arr, plot_m, plot_h, plot_j

def INa_func(parm):
    i = 0
    for Va in Vcm:
        mhj_sol = solve_ivp(INa_dt, tSpan, ini_con, args= (Va, parm), teval = tEval, method= Mthd)
#         Vmem = [Va if (S1_s <= k < S1_e) or (S2_s <= k < S2_e) or (S3_s <= k < S3_e) else Vhold for k in mhj_sol.t]
        Vmem = Vstim(mhj_sol.t, Va, frq) 
#         plt.plot(mhj_sol.t, Vmem)
        INa = parm[0]*(mhj_sol.y[0]**3)*mhj_sol.y[1]*mhj_sol.y[2] * (Vmem - Ena)
        plot_INa[i, :len(INa)] = INa
        time[i, :len(mhj_sol.t)] = mhj_sol.t
        size_arr[i] = np.uint16(len(mhj_sol.t))
        i = i+1
    for j in np.arange(len(Vcm)):
        INa_IVs[j] = np.min(plot_INa[j, : size_arr[j]])
    return INa_IVs


def INa_func_time(parm):
    i = 0
    for Va in Vcm:
        mhj_sol = solve_ivp(INa_dt, tSpan, ini_con, args= (Va, parm), teval = tEval, method= Mthd)
#         Vmem = [Va if (S1_s <= k < S1_e) or (S2_s <= k < S2_e) or (S3_s <= k < S3_e) else Vhold for k in mhj_sol.t]
        Vmem = Vstim(mhj_sol.t, Va, frq) 
#         plt.plot(mhj_sol.t, Vmem)
        INa = parm[0]*(mhj_sol.y[0]**3)*mhj_sol.y[1]*mhj_sol.y[2] * (Vmem - Ena)
        plot_INa[i, :len(INa)] = INa
        time[i, :len(mhj_sol.t)] = mhj_sol.t
        size_arr[i] = np.uint16(len(mhj_sol.t))
        i = i+1
    return time, plot_INa, size_arr
