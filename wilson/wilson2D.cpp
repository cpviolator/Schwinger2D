#include "schwinger2d_internal.h"
#include "utils.h"
//#include "measurements.h"
//#include "dirac_op.h"
//#include "inverters.h"
#include "hmc.h"

int main(int argc, char **argv) {
  
  param_t p;

  int Nx=32, Ny=32, Nd=2;

  p.Nx = Nx;
  p.Ny = Ny;
  p.beta = atof(argv[1]); 
  p.iter_hmc = atoi(argv[2]);
  p.therm = atoi(argv[3]);
  p.skip = atoi(argv[4]);
  p.chkpt = atoi(argv[5]);
  p.checkpoint_start = atoi(argv[6]);  
  p.n_step = atoi(argv[7]);
  p.tau = atof(argv[8]);
  
  p.smear_iter = atoi(argv[9]);
  p.alpha = atof(argv[10]);  
  p.seed = atoi(argv[11]);
  p.dynamic = (atoi(argv[12]) == 0 ? false : true);
  p.m = atof(argv[13]);
  p.max_iter_cg = atoi(argv[14]);
  p.eps = atof(argv[15]);
  
  //VOATOL params
  p.n_kr = 64;
  p.n_ev = 16;
  p.eig_tol = atof(argv[16]);
  p.eig_max_restarts = atoi(argv[17]);
  p.poly_acc = (atoi(argv[18]) == 0 ? false : true);
  p.amax = atof(argv[19]);
  p.amin = atof(argv[20]);
  p.poly_deg = atoi(argv[21]);

  //Measurements
  p.meas_pl = (atoi(argv[22]) == 0 ? false : true);
  p.meas_wl = (atoi(argv[23]) == 0 ? false : true);
  p.meas_pc = (atoi(argv[24]) == 0 ? false : true);
  p.meas_vt = (atoi(argv[25]) == 0 ? false : true);

  //Pseudo RNG seed
  srand48((long)p.seed);
  
  //Topology
  double top = 0.0;
  int top_int = 0;
  int top_old = 0;
  int top_stuck = 0;

  int histL = 101;
  std::vector<int> histQ(histL);
  double plaqSum = 0.0;
  int index = 0;
  for(int i = 0; i < histL; i++) histQ[i] = 0;
    
  int count = 0;
  string name;
  fstream outPutFile;

  int hmccount = 0;
  double expdHAve = 0.0;
  double dHAve = 0.0;

  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  printParams(p);
  
  field<Complex> *gauge = new field<Complex>(p);
  gaussStart(gauge);  // hot start

  //Start simulation
  double time0 = -((double)clock());
  int iter_offset = 0;
  int iter = 0;
  cout << setprecision(16);
  
  if(p.checkpoint_start > 0) {

    //Read in gauge field if requested
    //---------------------------------------------------------------------
    name = "gauge/gauge";
    constructName(name, p);
    name += "_traj" + to_string(p.checkpoint_start) + ".dat";	
    readGauge(gauge,name);
    iter_offset = p.checkpoint_start;    
  } else {
    
    //Thermalise from random start
    //---------------------------------------------------------------------
    for(iter=0; iter<2*p.therm; iter++){  
      //Perform HMC step
      accept = hmc(gauge, iter, expdHAve, dHAve, hmccount);
      
      double time = time0 + clock();
      cout << fixed << iter+1 << " ";              //Iteration
      cout << time/CLOCKS_PER_SEC << " " << endl;  //Time
    }    
    iter_offset = 2*p.therm;    
  }
  
  double beta0 = p.beta;
  bool measure = false;
  int freq = 200;
  int pulse = 50;
  int interval = 10;
  int strength = 2;
  
  //Begin thermalised trajectories
  //---------------------------------------------------------------------
  for(iter=iter_offset; iter<p.iter_hmc + iter_offset; iter++){
    
    //Perform HMC step
    accept = hmc(gauge, iter, expdHAve, dHAve, hmccount);    
    accepted += accept;
    
    //Measure the topological charge if trajectory is accepted
    //---------------------------------------------------------------------
    if(accept == 1) {
      /*
      top = measTopCharge(gauge, p);
      top_int = round(top);
      name = "data/top/top_charge";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %d\n", iter, top_int);
      fclose(fp);
      
      index = top_int + (histL-1)/2;
      histQ[index]++;
      if(top_old == top_int) top_stuck++;
      top_old = top_int;
      */
    }
    
    //Perform Measurements
    //---------------------------------------------------------------------
    if( (iter+1)%p.skip == 0) {
      
      count++; //Number of measurements taken

      //Checkpoint the gauge field?
      if( (iter+1)%p.chkpt == 0) {	  
	name = "gauge/gauge";
	constructName(name, p);
	name += "_traj" + to_string(iter+1) + ".dat";
	writeGauge(gauge, name);
      }
      
      //Plaquette action
      double plaq = measPlaq(gauge).real();
      plaqSum += plaq;

      //Dump simulation data to stdout
      double time = time0 + clock();
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << time/CLOCKS_PER_SEC << " ";                 //Time
      cout << plaqSum/count << " ";                       //Action
      cout << (double)top_stuck/(accepted) << " ";        //P(stuck)
      cout << expdHAve/hmccount << " ";                   //Average exp(-dH)
      cout << dHAve/hmccount << " ";                      //Average dH
      cout << (double)accepted/(count*p.skip) << " ";     //Acceptance
      cout << (double)p.beta << " ";                      //Current beta
      cout << top_int << endl;                            //T charge
	
      //Dump simulation data to file
      name = "data/data/data"; //I cannot make bricks without clay!
      constructName(name, p);
      name += ".dat";	
      sprintf(fname, "%s", name.c_str());	
      fp = fopen(fname, "a");	
      fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e %.16e %d\n",
	      iter+1,
	      time/CLOCKS_PER_SEC,
	      plaqSum/count,
	      (double)top_stuck/(accepted),
	      expdHAve/hmccount,
	      dHAve/hmccount,
	      (double)accepted/(count*p.skip),
	      top_int);
      fclose(fp);
      
      //Update topoligical charge histogram
      name = "data/top/top_hist";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "w");
      for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[i]);
      fclose(fp);

      //Physical observables
      //-------------------------------------------------------------      
      //Polyakov Loops      
      //if(p.meas_pl) measPolyakovLoops(gauge, iter, p);
      
      //Creutz Ratios
      //if(p.meas_wl) measWilsonLoops(gauge, plaq, iter, p);

      //Pion Correlation
      //if(p.meas_pc) measPionCorrelation(gauge, top_old, iter, p);

      //Vacuum Trace
      //if(p.meas_vt) measVacuumTrace(gauge, top_old, iter, p);
      //-------------------------------------------------------------
    }
  }
  return 0;
}
//-------------------------------------------------------------------------------
