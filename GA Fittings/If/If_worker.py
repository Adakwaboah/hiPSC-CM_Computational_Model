#If_worker function for the multiprocessing

import numpy as np

#constants
dt = 10 #time steps
tStart = 0
tEnd = 1500
nStep = np.int(np.ceil((tEnd - tStart) / dt))  # number of steps
t = np.linspace(tStart, tEnd, nStep + 1)
Vcm = np.arange(-110, -39, 10)
Vhold = -50  #mV

K_o = 5.4 #mM
Na_o = 140 #mM

K_i = 136.9896086978 #mM
Na_i = 15  #mM
R = 8.314 #J K^-1 mol^-1
T = 310 #K
F = 96.4867 #C.mmol^-1
Gfk = 0.0234346
Gfna = 0.0145654
Cmem = 1.38e-10
Ef = -17.0
Vm = np.zeros((nStep + 1, 1), dtype=np.float)
plot_y = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_tau_y = np.zeros((nStep, len(Vcm)), dtype=np.float)

plot_If = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_time = np.empty((nStep, 1), dtype=np.float)
plot_Im = np.empty((nStep, len(Vcm)), dtype=np.float)

def If_func(If_par):

    Ek = ((R*T)/F)*np.log(K_o/K_i)
    Ena = ((R*T)/F)*np.log(Na_o/Na_i)

    i = 0 #column counter

    for Vc in Vcm:
        tNow = tStart
        Vm[0] = Vhold
        # initial condition
        y = 0.0
        for iStep in np.arange(nStep):
            Vm[iStep+1] = (tNow < 100)* Vhold + (100 <= tNow < 600 )*Vc + (tNow >= 600 )*Vhold
            yss = 1/(1 + np.exp((Vm[iStep] + 80.6)/6.8))
            alpha_y = np.exp(-(If_par[2] + (If_par[3]*Vm[iStep])))
            beta_y = np.exp(If_par[4] + (If_par[5]*Vm[iStep]))
            tau_y = If_par[6]/(alpha_y + beta_y)
            dyf = ((yss - y)/ tau_y)*dt
            if_k = If_par[1] * y * (Vm[iStep] - Ek)
            if_na = If_par[0] * y * (Vm[iStep] - Ena)
            If = if_k + if_na            
            
            plot_y[iStep, i] = y
            plot_If[iStep, i] = If
            plot_tau_y[iStep, i] = tau_y
            plot_time[iStep] = tNow

            y = y+dyf
            tNow = tStart + iStep * dt
        i = i + 1

    If_peaks = np.min(plot_If, axis=0)
    return If_peaks

def If_func_time(If_par):
    Ek = ((R*T)/F)*np.log(K_o/K_i)
    Ena = ((R*T)/F)*np.log(Na_o/Na_i)

    i = 0 #column counter

    for Vc in Vcm:
        tNow = tStart
        Vm[0] = Vhold
        # initial condition
        y = 0.0
        for iStep in np.arange(nStep):
            Vm[iStep+1] = (tNow < 100)* Vhold + (100 <= tNow < 600 )*Vc + (tNow >= 600 )*Vhold
            yss = 1/(1 + np.exp((Vm[iStep] + 80.6)/6.8))
            alpha_y = np.exp(-(If_par[2] + (If_par[3]*Vm[iStep])))
            beta_y = np.exp(If_par[4] + (If_par[5]*Vm[iStep]))
            tau_y = If_par[6]/(alpha_y + beta_y)
            dyf = ((yss - y)/ tau_y)*dt

            if_k = If_par[1] * y * (Vm[iStep] - Ek)
            if_na = If_par[0] * y * (Vm[iStep] - Ena)
            If = if_k + if_na
            
            plot_y[iStep, i] = y
            plot_If[iStep, i] = If
            plot_tau_y[iStep, i] = tau_y
            plot_time[iStep] = tNow

            y = y+dyf
            tNow = tStart + iStep * dt
        i = i + 1
    If_peaks = np.min(plot_If, axis=0)
    return plot_time, plot_If, plot_y