#include "schwinger2d_internal.h"
#include "utils.h"
#include "hmc.h"
#include "iram.h"
#include "io.h"

int main(int argc, char **argv) {

  struct timeval start, end, total_start, total_end;
  gettimeofday(&total_start, NULL);
  double t_hmc = 0.0;
  double t_meas = 0.0;
  double t_total = 0.0;

  Param p;
  //Process Command line arguments
  for (int i=1; i<argc; i++){
    if(p.init(argc, argv, &i) == 0){
      continue;
    }
    printf("ERROR: Invalid option: %s\n", argv[i]);
    p.usage(argv);
    exit(0);
  }

  int k=0;
  if(argc > 1) p.init(argc, argv, &k);
    
  if(p.loop_max > std::min(p.Nx/2, p.Ny/2) && (p.meas_pl || p.meas_wl)) {
    cout << "Warning: requested Wilson loop max " << p.loop_max << " greater than ";
    cout << min(p.Nx/2, p.Ny/2) << ", truncating." << endl;
    p.loop_max = std::min(p.Nx/2, p.Ny/2);
  }
  
  if(p.sampler == S_MCHMC){
    double d = p.Nx * p.Ny * 2;
    cout << "Warning: manually changing tau from " << p.tau;
    p.tau = p.beta_eps * sqrt(d);
    p.n_step = 1;
    p.skip = 1;
      
    cout << " to " << p.tau << ", n_step from " << p.n_step << " to 1 and skip from " << p.skip << " to 1." << endl;
  }
  
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

  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  p.print();
  
  field<Complex> *gauge = new field<Complex>(p);
  gaussStart(gauge);  // hot start

  //Start simulation
  double time0 = -((double)clock());
  int iter_offset = 0;
  int iter = 0;
  cout << setprecision(16);
  
  HMC *HMCStep = new HMC(p);

  //Thermalise
  //---------------------------------------------------------------------
  gettimeofday(&start, NULL);  
  for(iter=p.checkpoint_start; iter<p.therm; iter++){      
      
    //Perform HMC step
    accept = HMCStep->hmc(gauge, iter);
    gettimeofday(&total_end, NULL);  
    t_total = ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
    cout << fixed << setprecision(16) << iter +1 << " "; //Iteration
    cout << t_total << " " << endl;                      //Time   
  }
  gettimeofday(&end, NULL);  
  t_hmc += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
  //Begin thermalised trajectories
  //---------------------------------------------------------------------
  for(;iter<p.iter_hmc + p.therm; iter++){
    
    //Measure the topological charge at each step if trajectory is accepted
    //---------------------------------------------------------------------
    top = measTopCharge(gauge);
    top_int = round(top);
    name = "data/top/top_charge";
    constructName(name, p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());
    fp = fopen(fname, iter==iter_offset ? "w" : "a");
    fprintf(fp, "%d %d\n", iter, top_int);
    fclose(fp);
    
    index = top_int + (histL-1)/2;
    histQ[index]++;
    if(top_old == top_int) top_stuck++;
    top_old = top_int;      

    //Plaquette action
    double plaq = HMCStep->measPlaq(gauge).real();
    plaqSum += plaq;
    
    gettimeofday(&total_end, NULL);
    t_total = ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
    
    //Dump simulation data to file
    name = "data/data/data"; //I cannot make bricks without clay!
    constructName(name, p);
    name += ".dat";	
    snprintf(fname, 100, "%s", name.c_str());	
    fp = fopen(fname, "a");	
    fprintf(fp, "%.06d %.16e %.16e %+.16e %+.16e %d %d\n",
	    iter,
	    t_total,
	    plaq,
	    HMCStep->exp_dH,
	    HMCStep->dH,
	    accept,
	    top_int);
    fclose(fp);
      
    //Update topoligical charge histogram
    name = "data/top/top_hist";
    constructName(name, p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());
    fp = fopen(fname, "w");
    for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[i]);
    fclose(fp);    
    
    //Perform HMC step
    gettimeofday(&start, NULL);
    accept = HMCStep->hmc(gauge, iter);
    accepted += accept;
    gettimeofday(&end, NULL);  
    t_hmc += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

    // Do a reversibility check
    if((iter+1)%p.reverse == 0 && p.sampler == S_HMC) {	  
      
      field<Complex> *gauge_old = new field<Complex>(p);
      blas::copy(gauge_old, gauge);
      bool reverse = HMCStep->hmc_reversibility(gauge_old, iter);
      if(!reverse) {
	cout << "Error in reversibility" << endl;
	exit(0);
      }
      
      // gauge_old should now be the same as gauge
      blas::axpy(-1.0, gauge, gauge_old);
      cout << "L2 norm of gauge(t0) - evolve_gauge(t0->t1->t0) = " << std::scientific << blas::norm2(gauge_old)/(p.Nx * p.Ny * 2) << endl;
      
      blas::copy(gauge_old, gauge);	
      double wilson_flow_tau = 1.0;
      int wilson_flow_steps = 50;
      double dt = wilson_flow_tau/wilson_flow_steps;
      for(int i=0; i<wilson_flow_steps; i++) wilsonFlow(gauge_old, dt);
      for(int i=0; i<wilson_flow_steps; i++) wilsonFlow(gauge_old, -dt);

      // gauge_old should now be the same as gauge
      blas::axpy(-1.0, gauge, gauge_old);
      cout << "L2 norm of gauge(t0) - flowed_gauge(t0->t1->t0) = " << std::scientific << blas::norm2(gauge_old)/(p.Nx * p.Ny * 2) << endl;
    }
    
    //Perform Measurements
    //---------------------------------------------------------------------
    if((iter+1)%p.skip == 0 && iter >= p.therm) {
      
      count++; //Number of measurements taken

      //Dump simulation data to stdout
      gettimeofday(&total_end, NULL);  
      cout << fixed << setprecision(16) << iter +1 << " ";// Iteration
      cout << t_total << " ";                             // Time
      cout << plaq << " ";                                // Action
      cout << plaqSum/(count*p.skip) << " ";              // Average Action
      cout << HMCStep->exp_dH << " ";                     // exp(-dH)
      cout << HMCStep->exp_dH_ave/(count*p.skip) << " ";  // Average exp(-dH)
      printf("%+.16f ", HMCStep->dH);                     // dH
      cout << HMCStep->dH_ave/(count*p.skip) << " ";      // Average dH
      cout << (double)accepted/(count*p.skip) << " ";     // Average Acceptance
      cout << top_int << endl;                            // T charge
      
      //Physical observables
      //-------------------------------------------------------------      
      //Polyakov Loops      
      //if(p.meas_pl) measPolyakovLoops(gauge, iter, p);
      
      //Creutz Ratios (for string tension)
      if(p.meas_wl) measWilsonLoops(gauge, plaq, iter);

      //Pion Correlation
      if(p.meas_pc) measPionCorrelation(gauge, iter);

      //Vacuum Trace
      //if(p.meas_vt) measVacuumTrace(gauge, top_old, iter, p);
      //-------------------------------------------------------------
    }    
  }
  return 0;
}
//-------------------------------------------------------------------------------
