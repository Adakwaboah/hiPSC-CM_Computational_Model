//----------------------------------------------------------------------------
// Akwaboah et al. Front. Physiol. 2021. doi: 10.3389/fphys.2021.675867 
// This file contains the ionic current formulations and initial conditions
// Called from the hipsc_singlecell_main.cc file
// Last modified 05/25/2021 - Makarand Deo
//----------------------------------------------------------------------------

#include <iostream> 
#include <iomanip> 
#include <math.h> 
#include <fstream> 
#include <cstdlib> 
#include <stdio.h>
#include "hipsc_const.h"

#define InitConditions "hipsc_steady_states.d"

using namespace std; 

class state_variables
{	
public:
 
  double Vm;
  double  m;
  double  h;
  double  j;
  double  Na_i;
  double  K_i;
  double  xtof;
  double  ytof;
  double  xtos;
  double  ytos;
  double  paF;
  double  paS;
  double pi;
  double  yf;
  double  d;
  double  f;
  double  fca;
  double f_tc;
  double f_tmc;
  double f_tmm;
  double f_cmi;
  double f_cms;
  double f_cq;
  double Ca_i;
  double Ca_sub;
  double Ca_up;
  double Ca_rel;
  double xs;
  double ua;
  double ui;
  
  
 double icallog;
 double icablog;
 double ipcalog;
 double incxlog;
 double inalog;
 double inaklog;
 double iktolog;
 double ik1log;
 double ikslog;
 double ikurlog;
 double ikrlog;
 double jsercalog;
 double jrellog;
 double inablog;
 double ikachlog;
 double iflog;
 double jcadifflog;
 double jtrlog;
    
 
};

const int count_init = 30;
double state[30];
double Iion;
double INa, Ito, Itof, Itos, IKr, If, INaK, ICaL, INaCa, IbCa, IbNa, IK1, IKs, IKur, IKAch, IpCa;
double ENa, EK, ECaL;
double am, bm, ah, bh, aj, bj, taum, tauh, tauj, mss, hss, jss;
double xtoss, ytoss, tau_xf, tau_yf, tau_xs, tau_ys;
double pa_ss, pi_ss, tau_paF, tau_paS, tau_pi, gkr;
double yss, alpha_y, beta_y, tau_y, If_k, If_na;
double d_ss, f_ss, fca_ss, alpha_dl, beta_dl, tau_d, tau_f, tau_fca;
double alpha_xs, beta_xs, tau_xs1, xs_ss, EKs;
double gkur, alpha_ua, beta_ua, tau_ua, ua_ss, alpha_ui, beta_ui, tau_ui, ui_ss;
double sigma, f_NaK;
double INaCa_num, INaCa_den;
double aK1, bK1, K1_ss;
double Irel=0, Fn, tau_u, u_ss, tau_v, v_ss, tau_w, w_ss,  Itr, Iup, Iupleak, caiont;
double Ca_cmdn, Ca_trpn, Ca_csqn, B1, B2;	
double dftc_dt, dftmc_dt, dftmm_dt, dfcmi_dt, dfcms_dt, dfcq_dt;
double jCa_diff, jrel, jup, jtr;



/************************************************************************/
/*  -------- function prototypes ----------------- */
void Get_state_variables_initial_condition();
void Assign_node_initial_condition(state_variables &rN);
double Total_transmembrane_currents(state_variables &rN, double &Ist, int location);

/*************************************************************************/
//void Get_state_variables_initial_condition();
void Assign_node_initial_condition(state_variables &rN);
double Total_transmembrane_currents(state_variables &rN, double &Ist, int location);

/*************************************************************************/
/* Read initial conditions (total 57) from file - Deo_mouse_pkj_steady_states.d  */

void Get_state_variables_initial_condition()
{
	FILE* Initptr;        
	Initptr = fopen(InitConditions, "r");	
	if (Initptr==NULL)
    {
       	cerr<<"cant open the file"<<endl;    	
   	}

	for(int jj=0; jj<count_init; jj++)
	{
	  	fscanf(Initptr,"%le",&state[jj]);
	  	//fscanf(Initptr,"%le",&ini_con[jj]);
	}
	
 	fclose(Initptr);
}

/**************************************************************************/
/*  Assign initial conditions to all the variables   */

void Assign_node_initial_condition(state_variables &rN)
{
    rN.Vm = state[0];
    rN.m = state[1];
    rN.h = state[2];
    rN.j = state[3];
    rN.Na_i = state[4];
    rN.K_i = state[5];
    rN.xtof = state[6];
    rN.ytof = state[7];
    rN.xtos = state[8];
    rN.ytos = state[9];
    rN.paF = state[10];
    rN.paS = state[11];
    rN.pi = state[12];
    rN.yf = state[13];
    rN.d = state[14];
    rN.f = state[15];
    rN.fca = state[16];
    rN.f_tc = state[17];
    rN.f_tmc = state[18];
    rN.f_tmm = state[19];
    rN.f_cmi = state[20];
    rN.f_cms = state[21];
    rN.f_cq = state[22];
    rN.Ca_i = state[23];
    rN.Ca_sub = state[24];
    rN.Ca_up = state[25];
    rN.Ca_rel = state[26];
    rN.xs = state[27];
    rN.ua = state[28];
    rN.ui = state[29];
           
}	
         

/**************************************************************************/
/*   Main function to calculate all ionic currents    */
double Total_transmembrane_currents(state_variables &rN, double &Istim)

{

	ENa = (R * T / (F)) * log(Na_o / rN.Na_i);
	EK = (R * T / (F)) * log(K_o / rN.K_i);
	ECaL = 0.6*(R * T / (2*F)) * log(Ca_o / rN.Ca_i);


	// Fast sodium current ____________________________________________________
	am = INa_par[1]*(rN.Vm+47.13)/(1-exp(-0.1*(rN.Vm+47.13)));
    bm = INa_par[2]*exp(-rN.Vm/11);
    
    if (rN.Vm < -40){
	
        ah = INa_par[3]*exp((80+rN.Vm)/-6.8);
        bh = INa_par[4]*exp(0.079*(rN.Vm))+INa_par[5]*exp(0.35*(rN.Vm));
        aj = (-INa_par[6]*exp(INa_par[7]*rN.Vm)-INa_par[8]*exp(-INa_par[9]*rN.Vm))*((rN.Vm+37.78)/(1+exp(INa_par[10]*(rN.Vm+79.23))));
        bj = (INa_par[11]*exp(-INa_par[12]*rN.Vm))/(1+exp(-0.1378*(rN.Vm+40.14)));
    }else{
	    ah = 0;
        bh = 1/(INa_par[13]*(1+exp((rN.Vm+INa_par[14])/-INa_par[15])));
        aj = 0;
        bj = (INa_par[16]*exp(-INa_par[17]*rN.Vm))/(1+exp(-INa_par[18]*(rN.Vm+INa_par[19])));
    }
    taum = 1/(am+bm);
    tauh = 1/(ah+bh);
    tauj = 1/(aj+bj);
    mss = am*taum;
    hss = ah*tauh;
    jss = aj*tauj;

    
    INa = K_INa*(INa_par[0] * (rN.m*rN.m*rN.m) * rN.h * rN.j * (rN.Vm - ENa));
   

	// Ito - Grandi-Pandit Formulation (GA Optimized) _____________________________________________
    xtoss = 1 / (1 + exp(-(rN.Vm - Ito_par[2]) / Ito_par[3]));
    ytoss = 1 / (1 + exp((rN.Vm + 41.1) / 6.68));
    //#fast Ito
    tau_xf = Ito_par[4]*(exp(-pow(((rN.Vm+Ito_par[5])/Ito_par[6]),2))) + Ito_par[7];
    tau_yf = Ito_par[8]*(exp(-pow((rN.Vm+Ito_par[9]),2)/Ito_par[10])) + Ito_par[11];
    Itof = Ito_par[0] * rN.xtof * rN.ytof * (rN.Vm - EK);
    
    //#slow Ito
    tau_xs = Ito_par[12]/(1+exp((rN.Vm+Ito_par[13])/Ito_par[14])) + Ito_par[15];
    tau_ys = Ito_par[16]/(1+exp((rN.Vm+Ito_par[17])/Ito_par[18])) + Ito_par[19];
    Itos = Ito_par[1] * rN.xtos * rN.ytos * (rN.Vm - EK);
    
	Ito = K_Ito*(Itof + Itos) ; //#25% allowable scaling adjustment

 
 	// IKr Kurata Formulation (GA Optimized)_________________________________________________
    pa_ss = 1.0/(1.0 + exp(-(rN.Vm + IKr_par[1])/IKr_par[2]));
    pi_ss = 1.0/(1.0 + exp((rN.Vm + IKr_par[3])/IKr_par[4]));
    tau_paF = IKr_par[5]/(IKr_par[6]*exp(rN.Vm/IKr_par[7]) + IKr_par[8]*exp(-rN.Vm/IKr_par[9]));
    tau_paS = IKr_par[5]/(IKr_par[10]*exp(rN.Vm/IKr_par[11]) + IKr_par[12]*exp(-rN.Vm/IKr_par[13]));
    tau_pi = 1/(IKr_par[14]*exp(-rN.Vm/IKr_par[15]) + IKr_par[16]*exp(rN.Vm/IKr_par[17]));
    gkr = IKr_par[0]*pow(K_o, 0.59);
    IKr = K_IKr*gkr*(rN.Vm - EK)*(0.6*rN.paF + 0.4*rN.paS)*rN.pi;


	//If Stewart Formulation (GA Optimized)__________________________________________________
    yss = 1/(1 + exp((rN.Vm + 80.6)/6.8));
    alpha_y = exp(-(If_par[2] + (If_par[3]*rN.Vm)));
    beta_y = exp(If_par[4] + (If_par[5]*rN.Vm));
    tau_y = If_par[6]/(alpha_y + beta_y);
    
    If_k = K_If*If_par[1] * rN.yf * (rN.Vm - EK);
    If_na = K_If*If_par[0] * rN.yf * (rN.Vm - ENa);
    If = If_k + If_na;
    if (If >0){
    	If=0;
	}
    
	// ICaL Kurata Formulation (GA Fitted)____________________________________________________
    d_ss = 1 / (1 + exp(-(rN.Vm + ICaL_par[1]) / ICaL_par[2])) ;
    f_ss = 1 / (1 + exp((rN.Vm + ICaL_par[3]) / ICaL_par[4])) ;
    fca_ss = K_mfca / (K_mfca + rN.Ca_sub);
    alpha_dl = -0.02839 * (rN.Vm + 35) / (exp(-(rN.Vm + 35) / 2.5) - 1) - (0.0849 * rN.Vm) / (exp(-(rN.Vm) / 4.8) - 1.00000001);
    beta_dl = 0.01143 * (rN.Vm - 5) / (exp((rN.Vm - 5) / 2.5) - 1);
    tau_d = 1 / (alpha_dl + beta_dl);
    tau_f = 257.1 * exp(-pow(((rN.Vm + 32.5) / 13.9),2)) + 44.3;
    tau_fca = fca_ss / alpha_fca;
    //ICaL = K_ICaL*ICaL_par[0] * (rN.Vm - ECaL) * rN.d * rN.f * rN.fca; //#25% allowable scaling adjustment
    ICaL = K_ICaL*ICaL_par[0]*rN.d*rN.f*rN.fca*4*(((rN.Vm-ICaL_par[5])*F)/Et)*(rN.Ca_i*exp((2*(rN.Vm-ICaL_par[5]))/Et) - 0.341*Ca_o)/(exp((2*(rN.Vm-ICaL_par[5]))/Et) - 1);

    
    //IKs Courtemanche Formulation______________________________________________________
    alpha_xs = (4e-5) * (rN.Vm - 19.9) / (1 - exp(-(rN.Vm - 19.9)/17));
    beta_xs = (3.5e-5) * (rN.Vm - 19.9) / (exp((rN.Vm - 19.9)/9) - 1);
    tau_xs1 = 1 / (2 * (alpha_xs + beta_xs));
    xs_ss = 1 / sqrt(1 + exp(-(rN.Vm - 19.9) / 12.7));
    EKs = ((R * T) / F) * log((0.01833 * Na_o + K_o) / (0.01833 * rN.Na_i + rN.K_i));
    IKs = K_IKs*G_ks * pow(rN.xs, 2) * (rN.Vm - EKs);
    
    
	// IKur Courtemanche Formulation___________________________________________________
    gkur = IKur_prm[0] + (IKur_prm[1] / (1 + exp(- (rN.Vm - IKur_prm[2]) / IKur_prm[3])));
    alpha_ua = IKur_prm[4] / (exp(-(rN.Vm + IKur_prm[5]) / IKur_prm[6]) + exp((-(rN.Vm - IKur_prm[7]) / IKur_prm[8])));
    beta_ua = IKur_prm[9] / (IKur_prm[10] + exp((rN.Vm + IKur_prm[11]) / IKur_prm[12]));
    tau_ua = 1 / (K_q10 * (alpha_ua + beta_ua));
    ua_ss = 1 / (1 + exp(-(rN.Vm + IKur_prm[13]) / IKur_prm[14]));
    alpha_ui = 1 / (IKur_prm[15] + exp(-(rN.Vm - IKur_prm[16]) / IKur_prm[17]));
    beta_ui = IKur_prm[18] * exp((rN.Vm - IKur_prm[19]) / IKur_prm[20]);
    tau_ui = 1 / (K_q10 * (alpha_ui + beta_ui));
    ui_ss = 1 / (1 + exp((rN.Vm - IKur_prm[21]) / IKur_prm[22]));
    IKur = K_IKur*gkur * pow(rN.ua,3) * rN.ui * (rN.Vm - EK);
     
    
	// INaK Courtemanche Formulation_____________________________________________________
    sigma = (1/7)*(exp(Na_o/67.3) - 1);
    f_NaK = 1/(1 + 0.1245*exp(-(0.1*(F*rN.Vm))/(R*T)) + 0.0365*sigma*exp(-(F*rN.Vm)/(R*T)));
    INaK = K_INaK*INaK_max * f_NaK * (1/(1+pow(K_mNai/rN.Na_i, 1.5))) * (K_o/(K_o + K_mKo));
    
    
	//INaCa Courtemanche Formulation______________________________________________________
    INaCa_num = INaCa_max * (exp(gamma_NaCa * (F / (R * T)) * rN.Vm) * pow(rN.Na_i,3) * Ca_o - exp((gamma_NaCa - 1) * (F / (R * T)) * rN.Vm) * pow(Na_o,3) * rN.Ca_i);
    INaCa_den = (pow(K_mNa,3) + pow(Na_o,3)) * (K_mCa + Ca_o) * (1 + k_sat * (exp((gamma_NaCa - 1) * (F / (R * T)) * rN.Vm)));
    INaCa = K_INaCa*(INaCa_num / INaCa_den); //#37
    
    
	//IKAch Kurata Formulation_____________________________________________________________
    IKAch = K_IKAch*gKAch*(rN.K_i - K_o*exp(-rN.Vm*F/(R*T)));


	//IK1 Grandi-Pandit Formulation_______________________________________________________________
    aK1 = 1.02 / (1 + exp(0.2385 * (rN.Vm - EK - 59.215)));
    bK1 = (0.49124 * exp(0.08032 * (rN.Vm + 5.476 - EK)) + exp(0.06175 * (rN.Vm - 594.31 - EK))) / (1 + exp(-0.5143 * (rN.Vm + 4.753 - EK)));
    K1_ss = aK1 / (aK1 + bK1);
    IK1 = K_IK1*0.35 * (sqrt(K_o / 5.4)) * K1_ss * (rN.Vm - EK); //#35

    //Background currents - Courtemanche___________________________________________________________
    IbCa = K_IbCa*G_bCa * (rN.Vm - ECaL);
    IbNa = K_IbNa*G_bNa * (rN.Vm - ENa);
   
    //#Ca2+ Pump Pump - maintains Ca_i to physiological levels
    IpCa = K_IpCa*IpCa_max * rN.Ca_i/(0.0005 + rN.Ca_i);
    
		
	//#Ca2+ ion homeostasis
    //#Ca Diffusion flux: Extracellular subspace and the intracellular
    jCa_diff = (rN.Ca_sub - rN.Ca_i)/tau_difCa;
   
           
    //#Ca handling by the SR - uptake, transfer and release
    jrel = Prel*(rN.Ca_rel - rN.Ca_sub)/(1 + pow((Krel/rN.Ca_sub),2)); //#SR release from JSR
    jup = Pup/(1 + Kup/rN.Ca_i);// #SR uptake by NSR
    jtr = (rN.Ca_up - rN.Ca_rel)/tau_tr; // # transfer from NSR to JSR
    
    	
	    
    //Sum of Ionic current density (already capacitance normalized)_________________________________
    Iion = INa + Ito + IKr + If + INaK + ICaL + INaCa + IbCa + IbNa + IK1 + IKs + IKur + IKAch + IpCa + Istim;
    



/*********************************************************************/
/* Ordinary Differential equations  */

// More efficient update of state variables which allows bigger time steps
rN.m = mss-(mss-rN.m)*exp(-dt/taum);
rN.h = hss-(hss-rN.h)*exp(-dt/tauh);
rN.j = jss-(jss-rN.j)*exp(-dt/tauj);


rN.xtof = xtoss-(xtoss-rN.xtof)*exp(-dt/tau_xf);
rN.ytof = ytoss-(ytoss-rN.ytof)*exp(-dt/tau_yf);
rN.xtos = xtoss-(xtoss-rN.xtos)*exp(-dt/tau_xs);
rN.ytos = ytoss-(ytoss-rN.ytos)*exp(-dt/tau_ys);


rN.paF = pa_ss-(pa_ss-rN.paF)*exp(-dt/tau_paF);
rN.paS = pa_ss-(pa_ss-rN.paS)*exp(-dt/tau_paS);
rN.pi = pi_ss-(pi_ss-rN.pi)*exp(-dt/tau_pi);

rN.yf = yss-(yss-rN.yf)*exp(-dt/tau_y);

rN.d   = d_ss - (d_ss-rN.d)*exp(-dt/tau_d);
rN.f   = f_ss - (f_ss-rN.f)*exp(-dt/tau_f);
rN.fca   = fca_ss - (fca_ss-rN.fca)*exp(-dt/tau_fca);

rN.xs = xs_ss-(xs_ss-rN.xs)*exp(-dt/tau_xs1);

rN.ua = ua_ss-(ua_ss-rN.ua)*exp(-dt/tau_ua);
rN.ui = ui_ss-(ui_ss-rN.ui)*exp(-dt/tau_ui);

// traditional update process
//rN.xtof += dt*((xtoss - rN.xtof)/ tau_xf);
//rN.ytof += dt*((ytoss - rN.ytof)/ tau_yf);
//rN.xtos += dt*((xtoss - rN.xtos)/ tau_xs);
//rN.ytos += dt*((ytoss - rN.ytos)/ tau_ys);

//rN.paF += dt*((pa_ss - rN.paF)/tau_paF);
//rN.paS += dt*((pa_ss - rN.paS)/tau_paS);
//rN.pi  += dt*((pi_ss - rN.pi)/tau_pi);

//rN.yf += dt*((yss - rN.yf)/ tau_y);

//rN.d   += dt*((d_ss-rN.d)/tau_d);
//rN.f   += dt*((f_ss-rN.f)/tau_f);
//rN.fca += dt*((fca_ss-rN.fca)/tau_fca);

//rN.xs  += dt*((xs_ss - rN.xs) / tau_xs1);

//rN.ua  += dt*((ua_ss - rN.ua) / tau_ua);
//rN.ui  += dt*((ui_ss - rN.ui) / tau_ui);


//#Ca buffering__________________________________
dftc_dt = (k_ftc*rN.Ca_i*(1 - rN.f_tc) - k_btc*rN.f_tc);
dftmc_dt = (k_ftmc*rN.Ca_i*(1 - rN.f_tmc - rN.f_tmm) - k_btmc*rN.f_tmc);
dftmm_dt = (k_ftmm*Mg_i*(1 - rN.f_tmc - rN.f_tmm) - k_btmm*rN.f_tmm);
dfcmi_dt = (k_fcm*rN.Ca_i*(1 - rN.f_cmi) - k_bcm*rN.f_cmi);
dfcms_dt = (k_fcm*rN.Ca_sub*(1 - rN.f_cms) - k_bcm*rN.f_cms);
dfcq_dt = (k_fcq*rN.Ca_rel*(1 - rN.f_cq) - k_bcq*rN.f_cq);

 //#intracellular ion concentrations___________________________________
rN.Ca_i += Cm*dt*((jCa_diff*V_sub - jup*V_up)/Vi - (CM_tot*dfcmi_dt + TC_tot*dftc_dt + TMC_tot*dftmc_dt));
rN.Ca_sub += Cm*dt*(jrel*(Vrel/V_sub) - (ICaL + IpCa + IbCa - 2*INaCa)/(2*F*V_sub) - (jCa_diff + CM_tot*dfcms_dt));
rN.Ca_rel += Cm*dt*(jtr - jrel - CQ_tot*dfcq_dt);
rN.Ca_up += Cm*dt*(jup - jtr*(Vrel/V_up));

rN.f_tc += dt*dftc_dt;
rN.f_tmc += dt*dftmc_dt;
rN.f_tmm += dt*dftmm_dt;
rN.f_cmi += dt*dfcmi_dt;
rN.f_cms += dt*dfcms_dt;
rN.f_cq += dt*dfcq_dt;

rN.Na_i  += dt*(-Cm*(3*INaK + 3*INaCa + INa + If_na + IbNa)/(F*Vi)) ;
rN.K_i   += dt*(-Cm*(Ito + IKr + If_k + IK1 + IKs + IKur + IKAch -2*INaK)/(F*Vi));
    
    
/* Write into current log variables	  */
rN.jsercalog = jup;
rN.jrellog = jrel;
rN.jcadifflog = jCa_diff;
rN.jtrlog = jtr; 

rN.icallog = ICaL;
rN.icablog = IbCa;
rN.ipcalog = IpCa;
rN.incxlog = INaCa;
rN.inalog  =  INa;
rN.inaklog = INaK;
rN.iktolog    = Ito;
rN.ik1log = IK1;
rN.ikslog = IKs;
rN.ikurlog = IKur;
rN.ikrlog = IKr;
rN.inablog = IbNa;
rN.ikachlog = IKAch;
rN.iflog = If;


return Iion;       /* returns total ionic current for each simulation which is further used for calculating transmembrane voltage in main file */
}
