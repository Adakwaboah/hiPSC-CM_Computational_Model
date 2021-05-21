#ICaL Worker function needed for the multiprocessing function

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

dt = 0.02#time steps
tSpan = (0, 980) #time start and end
tEval = np.arange(tSpan[0], tSpan[1], dt)

R = 8.3143 #J/(K.mol)
T = 310 #K
F = 96.4867 #C/mmol
Et = (R*T)/F
Ca_o = 2
# Ca_i_o = 0.019
d_o = 0
f_o = 1
fca_o = 1
ftc_o = 0.01
ftmc_o = 0.01
ftmm_o = 0.01
fcmi_o = 0.01
fcms_o = 0.01
fcq_o = 0.01
Ca_i_o = 1.02e-4
Ca_sub_o = 1.81e-3
Ca_up_o = 1.49
Ca_rel_o = 1.49

gCaL = 0.1238
K_ICaL = 1

#Ca_i constants
IpCa_max = 0.275 #pA/pF
Prel = 5*10.0 #ms-1 SR maximal release rate constant
Krel = 0.0012
Pup = 1.5*0.005 #0.005 #mM/ms maximal uptake rate constant
Kup = 0.0006#mM half-maximal Ca_up
tau_difCa = 0.04 #ms
tau_tr = 55*0.5 #ms

#Ca Buffering
Mg_i = 2.5 # intracellular Magnesium - used in buffering
k_ftc = 88.8
k_btmm = 0.751
k_btmc = 0.00751
k_bcq = 0.445
k_bcm = 0.542
k_fcm = 227.7
k_fcq = 0.534
k_ftmm = 2.277
k_ftmc = 227.7
k_btc = 0.446
TMC_tot = 0.062
TC_tot = 0.031
CM_tot = 0.045
CQ_tot = 10

Vcell = 20100 #um^3 converted from 3pl ie. 1000
V_sub = 0.01*Vcell
V_up = 0.0116*Vcell
Vi = 0.46*(Vcell - V_sub)
Vrel = 0.0012*Vcell
K_mfca = 0.00035
alpha_fca = 0.035

Vcm = np.arange(-40, 61, 10)
# Vcm = [-40,]
Vhold = -40
frq = 1

Mthd = 'BDF'

ini_con = (d_o, f_o, fca_o, ftc_o, ftmc_o, ftmm_o, fcmi_o, fcms_o, fcq_o, Ca_i_o, Ca_sub_o, Ca_up_o, Ca_rel_o, )

def Vm(t, Vcm, freq):
    T = 1000/freq
    if abs(t)%T >= T/2:
        return Vcm
    else:
        return Vhold
    
def Vstim(dur, Va, freq):
    Vstim = np.ones(len(dur))*Vhold
    T = 1000/freq
    i = 0
    for t in dur:
        if abs(t)%T >= T/2:
            Vstim[i] = Va
        i = i + 1
    return Vstim

def ICaL_dt(t, Y, Vc, ICaL_par):
    d = Y[0]
    f = Y[1]
    fca = Y[2]
    ftc = Y[3]
    f_tmc = Y[4]
    f_tmm = Y[5]
    f_cmi = Y[6]
    f_cms = Y[7]
    f_cq = Y[8]
    Ca_i = Y[9]
    Ca_sub = Y[10]
    Ca_up = Y[11]
    Ca_rel = Y[12]
    
    d_ss = 1/(1+ np.exp(-(Vm(t, Vc, frq) + ICaL_par[1])/ICaL_par[2]))
    f_ss = 1/(1+ np.exp((Vm(t, Vc, frq) + ICaL_par[3])/ICaL_par[4]))
    fca_ss = K_mfca/(K_mfca + Ca_i_o)
    alpha_d = -0.02839*(Vm(t, Vc, frq) + 35)/(np.exp(-(Vm(t, Vc, frq) + 35)/2.5) - 1) - (0.0849*Vm(t, Vc, frq))/(np.exp(-(Vm(t, Vc, frq))/4.8)-1.000000001)
    beta_d = 0.01143*(Vm(t, Vc, frq) - 5)/(np.exp((Vm(t, Vc, frq) - 5)/2.5) - 1)
    tau_d = 1/(alpha_d + beta_d)
    tau_f = 257.1*np.exp(-np.power((Vm(t, Vc, frq) + 32.5)/13.9, 2)) + 44.3
    tau_fca = fca_ss / alpha_fca
    ECaL = (R * T / (2*F)) * np.log(Ca_o / Ca_i)
    ICaL = ICaL_par[0]*d*f*fca*4*(((Vm(t, Vc, frq)-ICaL_par[5])*F)/Et)*(Ca_i*np.exp((2*(Vm(t, Vc, frq)-ICaL_par[5]))/Et) - 0.341*Ca_o)/(np.exp((2*(Vm(t, Vc, frq)-ICaL_par[5]))/Et) - 1)
    
    dd_dt = ((d_ss-d)/tau_d)
    df_dt = ((f_ss-f)/tau_f)
    dfca_dt = ((fca_ss-fca)/tau_fca)
    
    #Ca2+ Pump Pump - maintains Ca_i to physiological levels
    IpCa = IpCa_max * Ca_i/(0.0005 + Ca_i)
    
    #Ca2+ ion homeostasis
    #Ca Diffusion flux: Extracellular subspace and the intracellular
    jCa_diff = (Ca_sub - Ca_i)/tau_difCa
    
    #Ca handling by the SR - uptake, transfer and release
    jrel = Prel*(Ca_rel - Ca_sub)/(1 + (Krel/Ca_sub)**2) #SR release from JSR
    jup = Pup/(1 + Kup/Ca_i) #SR uptake by NSR
    jtr = (Ca_up - Ca_rel)/tau_tr # transfer from NSR to JSR
    
    #Ca buffering__________________________________
    dftc_dt = (k_ftc*Ca_i*(1 - ftc) - k_btc*ftc)
    dftmc_dt = (k_ftmc*Ca_i*(1 - f_tmc - f_tmm) - k_btmc*f_tmc)
    dftmm_dt = (k_ftmm*Mg_i*(1 - f_tmc - f_tmm) - k_btmm*f_tmm)
    dfcmi_dt = (k_fcm*Ca_i*(1 - f_cmi) - k_bcm*f_cmi)
    dfcms_dt = (k_fcm*Ca_sub*(1 - f_cms) - k_bcm*f_cms)
    dfcq_dt = (k_fcq*Ca_rel*(1 - f_cq) - k_bcq*f_cq)

    #intracellular ion concentrations___________________________________
    dCai_dt = ((jCa_diff*V_sub - jup*V_up)/Vi - (CM_tot*dfcmi_dt + TC_tot*dftc_dt + TMC_tot*dftmc_dt))
    dCasub_dt = (jrel*(Vrel/V_sub) - (ICaL)/(2*F*V_sub) - (jCa_diff + CM_tot*dfcms_dt))
    dCarel_dt = (jtr - jrel - CQ_tot*dfcq_dt)
    dCaup_dt = (jup - jtr*(Vrel/V_up))
    
    return dd_dt,df_dt,dfca_dt, dftc_dt, dftmc_dt, dftmm_dt, dfcmi_dt, dfcms_dt, dfcq_dt, dCai_dt, dCasub_dt, dCaup_dt, dCarel_dt

pts_limit = 30000
plot_ICaL = np.zeros((len(Vcm), pts_limit))
plot_Irel = np.zeros((len(Vcm), pts_limit))
plot_Fn = np.zeros((len(Vcm), pts_limit))
plot_d = np.zeros((len(Vcm), pts_limit))
plot_f = np.zeros((len(Vcm), pts_limit))
plot_fca = np.zeros((len(Vcm), pts_limit))
plot_Cai = np.zeros((len(Vcm), pts_limit))
plot_Caup = np.zeros((len(Vcm), pts_limit))
plot_Carel = np.zeros((len(Vcm), pts_limit))
plot_Casub = np.zeros((len(Vcm), pts_limit))
plot_jup = np.zeros((len(Vcm), pts_limit))
plot_jrel = np.zeros((len(Vcm), pts_limit))
plot_jtr = np.zeros((len(Vcm), pts_limit))
plot_jCa_diff = np.zeros((len(Vcm), pts_limit))

time = np.zeros((len(Vcm), pts_limit))
ICaL_IVs = np.zeros(len(Vcm))
size_arr = np.arange(len(Vcm))

def dffca(par):
    i = 0
    for Va in Vcm:
        dffca_sol = solve_ivp(ICaL_dt, tSpan, ini_con, args= (Va, par), teval = tEval, method= Mthd)
        Vmem = Vstim(dffca_sol.t, Va, frq)
        d = dffca_sol.y[0]
        f = dffca_sol.y[1]
        fca = dffca_sol.y[2]
        Ca_i = dffca_sol.y[9]
        Ca_sub = dffca_sol.y[10]
        Ca_up = dffca_sol.y[11]
        Ca_rel = dffca_sol.y[12]
        
        ECaL = (R * T / (2*F)) * np.log(Ca_o / Ca_i)
        ICaL = par[0]*d*f*fca*4*(((Vmem-par[5])*F)/Et)*(Ca_i - Ca_o*np.exp((-2*(Vmem-par[5]))/Et))/(1 - np.exp((-2*(Vmem-par[5]))/Et))
        jCa_diff = (Ca_sub - Ca_i)/tau_difCa
        jrel = Prel*(Ca_rel - Ca_sub)/(1 + (Krel/Ca_sub)**2) #SR release from JSR
        jup = Pup/(1 + Kup/Ca_i) #SR uptake by NSR
        jtr = (Ca_up - Ca_rel)/tau_tr # transfer from NSR to JSR
        
        plot_d[i, :len(dffca_sol.y[0])] = dffca_sol.y[0]
        plot_f[i, :len(dffca_sol.y[1])] = dffca_sol.y[1]
        plot_fca[i, :len(dffca_sol.y[2])] = dffca_sol.y[2]
        plot_ICaL[i, :len(dffca_sol.y[9])] = ICaL

        plot_Cai[i, :len(dffca_sol.y[9])] = Ca_i
        plot_Casub[i, :len(dffca_sol.y[10])] = Ca_sub
        plot_Caup[i, :len(dffca_sol.y[11])] = Ca_up
        plot_Carel[i, :len(dffca_sol.y[12])] = Ca_rel
        
        plot_jCa_diff[i, :len(dffca_sol.y[12])] = jCa_diff
        plot_jrel[i, :len(dffca_sol.y[12])] = jrel
        plot_jup[i, :len(dffca_sol.y[12])] = jup
        plot_jtr[i, :len(dffca_sol.y[12])] = jtr
        
        time[i, :len(dffca_sol.t)] = dffca_sol.t
        size_arr[i] = np.uint16(len(dffca_sol.t))
        i = i+1
    return time, size_arr, plot_d, plot_f, plot_fca, plot_Cai, plot_Casub, plot_Caup, plot_Carel, plot_jCa_diff, plot_jrel, plot_jup, plot_jtr, plot_ICaL

def ICaL_func(par):
    i = 0
    start_idx = 0
    for Va in Vcm:
        dffca_sol = solve_ivp(ICaL_dt, tSpan, ini_con, args= (Va, par), teval = tEval, method= Mthd)
        Vmem = Vstim(dffca_sol.t, Va, frq)
        d = dffca_sol.y[0]
        f = dffca_sol.y[1]
        fca = dffca_sol.y[2]
        Ca_i = dffca_sol.y[9]
        ECaL = (R * T / (2*F)) * np.log(Ca_o / Ca_i)
        ICaL = par[0]*d*f*fca*4*(((Vmem-par[5])*F)/Et)*(Ca_i*np.exp((2*(Vmem-par[5]))/Et) - 0.341*Ca_o)/(np.exp((2*(Vmem-par[5]))/Et) - 1)
        plot_ICaL[i, :len(ICaL)] = ICaL
        time[i, :len(dffca_sol.t)] = dffca_sol.t
        size_arr[i] = np.uint16(len(dffca_sol.t))
        i = i+1
    for j in np.arange(len(Vcm)):
        ICaL_IVs[j] = np.min(plot_ICaL[j, : size_arr[j]])
    return ICaL_IVs


def ICaL_func_time(par):
    i = 0
    for Va in Vcm:
        dffca_sol = solve_ivp(ICaL_dt, tSpan, ini_con, args= (Va, par), teval = tEval, method= Mthd)
        Vmem = Vstim(dffca_sol.t, Va, frq)
        d = dffca_sol.y[0]
        f = dffca_sol.y[1]
        fca = dffca_sol.y[2]
        Ca_i = dffca_sol.y[9]
        ECaL = (R * T / (2*F)) * np.log(Ca_o / dffca_sol.y[6])
        ICaL = par[0]*d*f*fca*4*(((Vmem-par[5])*F)/Et)*(Ca_i*np.exp((2*(Vmem-par[5]))/Et) - 0.341*Ca_o)/(np.exp((2*(Vmem-par[5]))/Et) - 1)
        plot_ICaL[i, :len(ICaL)] = ICaL
        time[i, :len(dffca_sol.t)] = dffca_sol.t
        size_arr[i] = np.uint16(len(dffca_sol.t))
        i = i+1
    return time, plot_ICaL, size_arr

