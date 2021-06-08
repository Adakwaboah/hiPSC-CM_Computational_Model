
//*************************************************************************//
//
// HiPSC-CM Single cell biophysical model
// Akwaboah et al. Front. Physiol. 2021. doi: 10.3389/fphys.2021.675867 
// 
// Required files: hipsc_ion.h
//                 hipsc_const.h
//				   hipsc_steady_states.d
//                
// Last modified 05/25/2021 - Makarand Deo
//*************************************************************************//

#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hipsc_ion.h"

//*************************************************************************//
/* User modified parameters  */


char IVfile  []  = "hipsc_IV.txt";        /* Output file with all currents and Vm */



/* Time parameters   */

double simulation_time = 60000.0;                  /* milliseconds - Total time of simulation  */
double frequency = 1;                            /* Hz - Frequency of stimulation  */
int pacing = 0;                                  // make 1 for enabling pacing; 0 for spontaneous firing 
int steadywrite = 0;                               /* make 1 to write steady states at the end of simulation.  */
      

char steadycond []  = "hipsc_steady_states.d";  /* Name of steady states file  */

/****************************************************************************/


/************** main function starts ****************************************/

int main(int argc, char *argv[])                                               
{

  state_variables * Node = new state_variables[1];
  
  int i, time_to_disp=0;
  FILE *fout,  *stateout;
  double Stim, sim_time;
  double Iion;
  double time_shift=0.0;
  double time_off=0.5;
  double delta_t=0.0;
  int set=1,counter=1;
  int done =0;
  int dtcnt = 0;
  
  if (pacing) {
	stim_mag=-50.0;
	K_If = 0.85*0.5;    // for disabling spontaneous activity 
	}
	else {
	stim_mag =0;
	K_If = 0.85;
	 }


  Get_state_variables_initial_condition();  /* Put initial conditions in the vector state[...] */
  Assign_node_initial_condition(Node[0]);
  
  fout = fopen(IVfile,"w");     
  
  cout << "Computation progress ..." << endl;
   
  /* Main simulation time loop starts    */
  for (sim_time=0.0; sim_time<simulation_time; sim_time+=dt)
  {
     /* Stimulus */
     //if ((sim_time >= time_shift + delta_t) && (sim_time <= time_off + time_shift +delta_t) && (sim_time <=10000))       /* use for stoping stimulus after certain time defined by last logical comparison  */
     if ((sim_time >= time_shift + delta_t) && (sim_time <= time_off + time_shift +delta_t))
     { 
        Stim = stim_mag;//-50; //0.0;
        
     } 
     else 
     {
          Stim = 0.0;
	      set =1;
     }
     
	 /*
	 // simulating Ach experiments: 0-20,000 ms control; 20,000-40,000 Ach exposure; >40,000 washout  
	 if ((sim_time >= 20000)&&(sim_time<40000))
     {
     	K_IKAch = 3;
	 }
	 else
	 {
	 	K_IKAch = 1;
	 }
	 */

     /* Main call to the ionic model    */
     Iion = Total_transmembrane_currents(Node[0], Stim); 
     Node[0].Vm += dt*-1.0*Iion;  /* Update Vm   */


     i = (simulation_time-sim_time)/dt;              /* integer to control time step for file writing  */
       
     /* Write variables to the files */
     if ((i % 100 == 0)&&(sim_time >0))
     {   
         
         fprintf(fout,"%8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf %8.12lf\n", sim_time, Node[0].Vm, Node[0].icallog, Node[0].ipcalog, Node[0].incxlog, Node[0].inalog, Node[0].inaklog, Node[0].iktolog, Node[0].ik1log, Node[0].ikslog, Node[0].ikurlog, Node[0].ikrlog, Node[0].jsercalog, Node[0].jrellog, Node[0].K_i, Node[0].Na_i, Node[0].Ca_i, Node[0].ikachlog, Node[0].iflog, Node[0].icablog, Node[0].inablog, Node[0].jcadifflog, Node[0].jtrlog, Node[0].Ca_up, Node[0].Ca_rel, Node[0].Ca_sub);
         
         
             
       }
       
     if (sim_time >= time_to_disp)
     {
        cout << time_to_disp << " ms" << endl;        /* displays time progress of simulation */
        time_to_disp+=10;
     }
      		    
     if((sim_time>=(1000.0*counter/frequency)+time_shift) && (set==1))
     {
        set=0;
        delta_t=0.0+(1000.0*counter/frequency);
        counter++;
     }


   
     /* save steady state values  */
     if(steadywrite)            /* saves the steady state file when steadywrite = 1,   */
     {
        
        if( (sim_time >= simulation_time-1000.0) && (done ==0) )
        {
            stateout = fopen(steadycond,"w");
        
            fprintf(stateout,"%8.12lf\n", Node[0].Vm);
            fprintf(stateout,"%8.12lf\n", Node[0].m);
            fprintf(stateout,"%8.12lf\n", Node[0].h);
            fprintf(stateout,"%8.12lf\n", Node[0].j);
            fprintf(stateout,"%8.12lf\n", Node[0].Na_i);
            fprintf(stateout,"%8.12lf\n", Node[0].K_i);
            fprintf(stateout,"%8.12lf\n", Node[0].xtof);
            fprintf(stateout,"%8.12lf\n", Node[0].ytof);
            fprintf(stateout,"%8.12lf\n", Node[0].xtos);
            fprintf(stateout,"%8.12lf\n", Node[0].ytos);
            fprintf(stateout,"%8.12lf\n", Node[0].paF);
            fprintf(stateout,"%8.12lf\n", Node[0].paS);
            fprintf(stateout,"%8.12lf\n", Node[0].pi);
            fprintf(stateout,"%8.12lf\n", Node[0].yf);
            fprintf(stateout,"%8.12lf\n", Node[0].d);
            fprintf(stateout,"%8.12lf\n", Node[0].f);
            fprintf(stateout,"%8.12lf\n", Node[0].fca);
            fprintf(stateout,"%8.12lf\n", Node[0].f_tc);
            fprintf(stateout,"%8.12lf\n", Node[0].f_tmc);
            fprintf(stateout,"%8.12lf\n", Node[0].f_tmm);
            fprintf(stateout,"%8.12lf\n", Node[0].f_cmi);
            fprintf(stateout,"%8.12lf\n", Node[0].f_cms);
            fprintf(stateout,"%8.12lf\n", Node[0].f_cq);
            fprintf(stateout,"%8.12lf\n", Node[0].Ca_i);
            fprintf(stateout,"%8.12lf\n", Node[0].Ca_sub);
            fprintf(stateout,"%8.12lf\n", Node[0].Ca_up);
            fprintf(stateout,"%8.12lf\n", Node[0].Ca_rel);
             fprintf(stateout,"%8.12lf\n", Node[0].xs);
            fprintf(stateout,"%8.12lf\n", Node[0].ua);
            fprintf(stateout,"%8.12lf\n", Node[0].ui);
           	
	            
            done++;
	        fclose(stateout);
         }
    }
          
          
  } /* main simulation time loop ends here */
    
  fclose(fout);

  cout << "IV file written: " << IVfile << endl;

  if (steadywrite)
  {
     cout << "Steady State file written: " << steadycond << endl;
     cout << " " << endl;
  }  

  cout << "DONE DONE DONE !!!" << endl;
 
  system("pause");
  return 0;

} /* End of main loop  */

