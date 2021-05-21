#Ito worker function for the multiprocessing function

import numpy as np

dt = 0.1  # time step
tStart = -15
tEnd = 80.0
nStep = np.int(np.ceil((tEnd - tStart) / dt))
Vcm = np.arange(-40, 51, 10)
Vh_pre = -80    # mV,prepulse Holding Voltage
Vh_brief = -50 #mV
Vh_post = -80   # mV post holding

Vm = np.zeros((nStep + 1, 1), dtype=np.float)
plot_Ito = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_xtof = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_ytof = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_xtos = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_ytos = np.zeros((nStep, len(Vcm)), dtype=np.float)
plot_time = np.empty((nStep, 1), dtype=np.float)

t = np.linspace(tStart, tEnd, nStep + 1)


def Ito_func(par):
    i = 0  #column selector

    for Vc in Vcm:
        tNow = tStart
        Vm[0] = Vh_pre
        # initial condition
        xtof = 0.0
        ytof = 1.0
        xtos = 0.0
        ytos = 1.0

        EK = -82.8  # mV, Potassium Nernst potential
        for iStep in np.arange(nStep):
            Vm[iStep+1] = (-10 < tNow < 0.0)* Vh_brief + (tNow <= -10)*Vh_pre + (0.0 <= tNow < 500.0) * Vc +  (tNow >= 500.0) * Vh_post
            xtoss = 1 / (1 + np.exp(-(Vm[iStep] - par[2]) / par[3]))
            ytoss = 1 / (1 + np.exp((Vm[iStep] + 41.1) / 6.68))
            
            #fast Ito
            tau_xf = par[4]*(np.exp(-((Vm[iStep]+par[5])/par[6])**2)) + par[7]
            tau_yf = par[8]*(np.exp(-((Vm[iStep]+par[9])**2)/par[10])) + par[11]
            dxtof = ((xtoss - xtof)/ tau_xf)*dt
            dytof = ((ytoss - ytof)/ tau_yf)*dt
            Itof = par[0] * xtof * ytof * (Vm[iStep] - EK)
            
            #slow Ito
            tau_xs = par[12]/(1+np.exp((Vm[iStep]+par[13])/par[14])) + par[15]
            tau_ys = par[16]/(1+np.exp((Vm[iStep]+par[17])/par[18])) + par[19]
            dxtos = ((xtoss - xtos)/ tau_xs)*dt
            dytos = ((ytoss - ytos)/ tau_ys)*dt
            Itos = par[1] * xtos * ytos * (Vm[iStep] - EK)
            
            Ito = Itof + Itos

            plot_Ito[iStep, i] = Ito
            plot_time[iStep] = tNow
            
            xtof = xtof+dxtof
            ytof = ytof+dytof
            xtos = xtos+dxtos
            ytos = ytos+dytos
            
            tNow = tStart + iStep * dt
        i = i + 1

    Ito_peaks = np.max(plot_Ito, axis=0)

    return Ito_peaks

def Ito_func_time(par):
    i = 0  #column selector

    for Vc in Vcm:
        tNow = tStart
        Vm[0] = Vh_pre
        # initial condition
        xtof = 0.0
        ytof = 1.0
        xtos = 0.0
        ytos = 1.0

        EK = -82.8  # mV, Potassium Nernst potential
        for iStep in np.arange(nStep):
            Vm[iStep+1] = (-10 < tNow < 0.0)* Vh_brief + (tNow <= -10)*Vh_pre + (0.0 <= tNow < 500.0) * Vc +  (tNow >= 500.0) * Vh_post
            xtoss = 1 / (1 + np.exp(-(Vm[iStep] - par[2]) / par[3]))
            ytoss = 1 / (1 + np.exp((Vm[iStep] + 41.1) / 6.68))
            
            #fast Ito
            tau_xf = par[4]*(np.exp(-((Vm[iStep]+par[5])/par[6])**2)) + par[7]
            tau_yf = par[8]*(np.exp(-((Vm[iStep]+par[9])**2)/par[10])) + par[11]
            dxtof = ((xtoss - xtof)/ tau_xf)*dt
            dytof = ((ytoss - ytof)/ tau_yf)*dt
            Itof = par[0] * xtof * ytof * (Vm[iStep] - EK)
            
            #slow Ito
            tau_xs = par[12]/(1+np.exp((Vm[iStep]+par[13])/par[14])) + par[15]
            tau_ys = par[16]/(1+np.exp((Vm[iStep]+par[17])/par[18])) + par[19]
            dxtos = ((xtoss - xtos)/ tau_xs)*dt
            dytos = ((ytoss - ytos)/ tau_ys)*dt
            Itos = par[1] * xtos * ytos * (Vm[iStep] - EK)
            
            Ito = Itof + Itos

            plot_Ito[iStep, i] = Ito
            plot_xtof[iStep, i] = xtof
            plot_xtos[iStep, i] = xtos
            plot_ytof[iStep, i] = ytof
            plot_ytos[iStep, i] = ytos
            
            plot_time[iStep] = tNow
            
            xtof = xtof+dxtof
            ytof = ytof+dytof
            xtos = xtos+dxtos
            ytos = ytos+dytos
            
            tNow = tStart + iStep * dt
        i = i + 1

    return plot_time, plot_Ito, plot_xtof, plot_ytof, plot_xtos, plot_ytos
