'''IKr worker file - Grandi-Pandit Model'''
import numpy as np

Vcm = np.arange(-40, 61, 20)    # mV, Command Voltage
Vh_pre = -80    # mV,prepulse Holding Voltage
Vh_post = -50   # mV post holding
dt = 0.2  # time step
splice_pt = np.uint16(100/dt)
tStart = -10
tEnd = 400.0
nStep = np.int(np.ceil((tEnd - tStart) / dt))  # number of steps

R = 8.314 #J K^-1 mol^-1
T = 310 #K
F = 96.4867 #C.mmol^-1
K_o = 5.4 #mM
K_i = 130 #mM
EK = ((R*T)/F)*np.log(K_o/K_i)
# Store the state variables

#     gKr = par[0] * np.sqrt(Ko/ 5.4)

Vm = np.zeros((nStep + 1, 1), dtype=np.float)
plot_paF = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_paS = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_pi = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_IKr = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_time = np.empty((nStep, 1), dtype=np.float)


t = np.linspace(tStart, tEnd, nStep + 1)

# [0.0294, 15, 22.4, 3e-4, 14.1, 5, 7.3898e-5, 3.3328, 5.1237, 14.1, 6.5]

def IKr_func(par):

    i = 0  #column selector

    for Vc in Vcm:
        tNow = tStart
        Vm[0] = Vh_pre
        paF = 0  # initial condition
        paS = 0
        pi = 1

        for iStep in np.arange(nStep):
            Vm[iStep+1] = (tNow < 0.0)*Vh_pre + (0.0 <= tNow < 300.0) * Vc +  (tNow >= 300.0) * Vh_post
            pa_ss = 1.0/(1.0 + np.exp(-(Vm[iStep] + par[1])/par[2]))
            pi_ss = 1.0/(1.0 + np.exp((Vm[iStep] + par[3])/par[4]))
            tau_paF = par[5]/(par[6]*np.exp(Vm[iStep]/par[7]) + par[8]*np.exp(-Vm[iStep]/par[9]))
            tau_paS = par[5]/(par[10]*np.exp(Vm[iStep]/par[11]) + par[12]*np.exp(-Vm[iStep]/par[13]))
            tau_pi = 1/(par[14]*np.exp(-Vm[iStep]/par[15]) + par[16]*np.exp(Vm[iStep]/par[17]))
            gkr = par[0]*np.power(K_o, 0.59)
            IKr = gkr*(Vm[iStep] - EK)*(0.6*paF + 0.4*paS)*pi
            
            dpaF = ((pa_ss - paF)/tau_paF)*dt
            dpaS = ((pa_ss - paS)/tau_paS)*dt
            dpi = ((pi_ss - pi)/tau_pi)*dt

            plot_paF[iStep, i] = paF
            plot_paS[iStep, i] = paS
            plot_pi[iStep, i] = pi
            plot_IKr[iStep, i] = IKr
            plot_time[iStep] = tNow
            paF = paF + dpaF
            paS = paF + dpaF
            pi = pi + dpi
            tNow = tStart + iStep * dt
        i = i + 1
    Ikr_peaks = np.max(plot_IKr[splice_pt:], axis=0)
    return Ikr_peaks

def IKr_func_time(par):

    i = 0  #column selector

    for Vc in Vcm:
        tNow = tStart
        Vm[0] = Vh_pre
        paF = 0  # initial condition
        paS = 0
        pi = 1
        for iStep in np.arange(nStep):
            Vm[iStep+1] = (tNow < 0.0)*Vh_pre + (0.0 <= tNow < 300.0) * Vc +  (tNow >= 300.0) * Vh_post
            pa_ss = 1.0/(1.0 + np.exp(-(Vm[iStep] + par[1])/par[2]))
            pi_ss = 1.0/(1.0 + np.exp((Vm[iStep] + par[3])/par[4]))
            tau_paF = par[5]/(par[6]*np.exp(Vm[iStep]/par[7]) + par[8]*np.exp(-Vm[iStep]/par[9]))
            tau_paS = par[5]/(par[10]*np.exp(Vm[iStep]/par[11]) + par[12]*np.exp(-Vm[iStep]/par[13]))
            tau_pi = 1/(par[14]*np.exp(-Vm[iStep]/par[15]) + par[16]*np.exp(Vm[iStep]/par[17]))
            gkr = par[0]*np.power(K_o, 0.59)
            IKr = gkr*(Vm[iStep] - EK)*(0.6*paF + 0.4*paS)*pi
            
            dpaF = ((pa_ss - paF)/tau_paF)*dt
            dpaS = ((pa_ss - paS)/tau_paS)*dt
            dpi = ((pi_ss - pi)/tau_pi)*dt

            plot_paF[iStep, i] = paF
            plot_paS[iStep, i] = paS
            plot_pi[iStep, i] = pi
            plot_IKr[iStep, i] = IKr
            plot_time[iStep] = tNow
            paF = paF + dpaF
            paS = paF + dpaF
            pi = pi + dpi
            tNow = tStart + iStep * dt
        i = i + 1
    Ikr_peaks = np.max(plot_IKr[splice_pt:], axis=0)
    return plot_time, plot_IKr, plot_paF, plot_paS, plot_pi