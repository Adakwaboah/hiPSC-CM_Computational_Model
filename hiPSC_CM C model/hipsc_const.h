//-------------------------------------------------------
//Akwaboah et al. Front. Physiol. 2021. doi: 10.3389/fphys.2021.675867 
// This file defines all the model constants
// Called from hipsc_ion.h file
// Last modified 05/25/2021 - Makarand Deo
//----------------------------------------------------

#include <iostream> 
#include <iomanip> 
#include <math.h> 
#include <fstream> 
#include <cstdlib> 
#include <stdio.h>


double dt = 0.0005;  // simulation step time

double R = 8.3143;// #J/(K.mol)
double T = 310; //#K
double F = 96.4867; // #C/mmol
double Cm = 100;//80; // #pF; average capacitance for INa
double Et = (R*T)/F;
//double V_i = 13668;// #um^3; intracellular volume

double K_If;   // Now initialized in the main .cc file
double stim_mag; // Initialized in the main .cc file 

//double K_If = 0.85*0.5; //0.85;   

//--------------------------------
// control model paramaters
double K_ICaL = 1.1; //#1.1  #5 ////////////////////////////
double K_INaK = 1.7*1.978*0.95;// #1.7*1.978  #1.7*1.79    #8
double K_INaCa = 1*9*1.8; // #1*9    #9
double K_IK1 = 0.18;// #0.18   #10   #########
double Pup = 0.005;//0.005;
//------------------------------------

/*
//--------------------------------
// Use these for isoproterenol simulations
double K_ICaL = 1.1*2.0; //#1.1  #5
double K_INaK = 1.7*1.978*0.95*1.3;// #1.7*1.978  #1.7*1.79    #8
double K_INaCa = 1*9*1.8*1.3; // #1*9    #9
double K_IK1 = 0.18*0.8;// #0.18   #10   #########
double Pup = 0.005*1.2;//0.005*1.2;// #1.5*0.005 #0.005 #mM/ms maximal uptake rate constant
//--------------------------------
*/
 
double Ca_o = 2;//2;// #mM    // Increase to 4 or 6 to simulate hypercalcemia
double K_o = 5.4;             // Increase to 10 or 12 to simulate hyperkalemia  

double K_INa = 1;//  #1
double K_Ito = 0.75*0.85;// #0.75  #2

double K_IKr = 1.25;// #3 

double K_IKur = 0.8;// #0.8     #6

double K_IKs = 0.3;//   #7
double K_IKAch = 1;//1;//    #11

double K_IbNa = 1.6;// #12
double K_IbCa = 1.6; // #13
double K_IpCa = 5.2;//    #14

//#ICaL and Ca_i constants
double IpCa_max = 0.275; //#pA/pF
double Prel = 5*2;// #ms-1 SR maximal release rate constant
double Krel = 0.0012;
//double Pup = 0.005;// #1.5*0.005 #0.005 #mM/ms maximal uptake rate constant
double Kup = 0.0006;//#mM half-maximal Ca_up
double tau_difCa = 0.04;// #ms
double tau_tr = 60;// #55*0.5 #ms # 60  #############

double prelold =5*2;


//#Ca Buffering
double Mg_i = 2.5;// # intracellular Magnesium - used in buffering
double k_ftc = 88.8;
double k_btmm = 0.751;
double k_btmc = 0.00751;
double k_bcq = 0.445;
double k_bcm = 0.542;
double k_fcm = 227.7;
double k_fcq = 0.534;
double k_ftmm = 2.277;
double k_ftmc = 227.7;
double k_btc = 0.446;
double TMC_tot = 0.062;
double TC_tot = 0.031;
double CM_tot = 0.045;
double CQ_tot = 10;
double K_mfca = 0.00035;
double alpha_fca = 0.035;

double Vcell = 3500; // #20100.0 #um^3 converted from 3pl ie. 1000#####
double V_sub = 0.01*Vcell;
double V_up = 0.0116*Vcell;
double Vi = 0.46*(Vcell - V_sub);
double Vrel = 0.0012*Vcell;

//#IKs
double G_ks = 0.129;// #nS/pF


//#IKur
double IKur_prm[23] = {0.005, 0.05, 15, 13, 0.65, 10, 8.5, 30, 59, 0.65, 2.5, 82, 17, 30.3, 9.6, 21, 185, 28, 1, 158, 16, 99.45, 27.48};
double K_q10 = 3;

double gKAch = 0.0011*(pow(K_o, 0.41));

//#INaK
double INaK_max = 4.3; //4.3; //#pA/pF
double K_mNai = 87.5; //#mM Nai half-saturation constant
double K_mKo = 1.0; //1.5; //#mM Ki half-saturation constant

//#INaCa
double gamma_NaCa = 0.35;
double k_sat = 0.1;
double K_mCa = 1.38; //#mM half-saturation constant
double K_mNa = 10; //#mM
double INaCa_max = 100; //10; //#pA/pF

//#Background currents
double G_bCa = 0.00113;
double G_bNa = 0.000674;
double G_bK = 0.000674;
double gK1 = 0.09; //#nS/pF

double Na_o = 140;// #mM   


double INa_par[20]= {1.01872382e+01, 2.81165049e-01, 4.11026646e-02, 1.85913230e-01,
       2.33559087e+00, 2.69343146e+05, 9.44262616e+04, 1.95789481e-01,
       4.34129546e-05, 5.95815316e-02, 4.29535950e-01, 1.26055107e-01,
       5.45387213e-03, 1.11372808e-01, 8.23823189e-01, 7.98343615e+00,
       4.37170307e-01, 2.87045755e-07, 2.42211297e-02, 5.22647290e+01};

double Ito_par[20]={1.07795129e-03, 1.55692740e-01, 1.46825076e+01, 1.59100592e+01,
       3.83118368e+00, 4.80328534e+01, 4.89599981e+01, 2.85524250e-01,
       1.23540497e+02, 5.03128630e+01, 3.16319213e+02, 1.32886316e+01,
       4.13205093e+00, 2.30054808e+00, 2.67131386e+01, 1.12598680e-01,
       1.13596188e+03, 4.95801337e+01, 5.53181883e+00, 1.97135761e+01};

double IKr_par[18]={1.86484701e-02, 6.42652882e+00, 1.18799783e+01, 4.01163188e+01,
       1.40822232e+01, 9.63756681e-01, 4.23722972e-02, 1.70356140e+01,
       1.24135620e-03, 2.33474212e+01, 4.56822754e-03, 1.46185713e+01,
       1.71859133e-04, 1.73714985e+01, 8.47065470e-02, 5.21505971e+01,
       7.59956413e-01, 1.07977148e+02};

double If_par[7]={5.02968027e-02, 7.02798870e-02, 1.68826991e+00, 3.67280874e-02,
       5.17628426e+00, 5.35050858e-02, 4.01031320e+03};

double ICaL_par[6]= { 0.18722844, 1.00592002, 6.38249172, 29.16060561,  3.60776178,  0.09457633};


