#include "schwinger2d_internal.h"
#include "utils.h"
#include "hmc.h"
#include "iram.h"
#include "io.h"

int main(int argc, char **argv) {

  struct timeval start, end, total_start, total_end;
  gettimeofday(&total_start, NULL);
  double t_hmc =  0.0;
  double t_meas = 0;
  double t_total =  0.0;
  
  param_t p;

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
  
  //eigensolver params
  p.deflate = (atoi(argv[16]) == 0 ? false : true);
  p.n_kr = atoi(argv[17]);
  p.n_ev = atoi(argv[18]);
  p.n_conv = atoi(argv[19]);
  p.eig_tol = atof(argv[20]);
  p.eig_max_restarts = atoi(argv[21]);
  p.poly_acc = (atoi(argv[22]) == 0 ? false : true);
  p.amax = atof(argv[23]);
  p.amin = atof(argv[24]);
  p.poly_deg = atoi(argv[25]);
  p.block_scheme[0] = atoi(argv[26]);
  p.block_scheme[1] = atoi(argv[27]);
  p.n_low = atoi(argv[28]);
  p.n_deflate = atoi(argv[29]); 
  
  //Measurements
  //p.meas_pl = (atoi(argv[22]) == 0 ? false : true);
  p.meas_wl = (atoi(argv[30]) == 0 ? false : true);
  p.meas_pc = (atoi(argv[31]) == 0 ? false : true);
  //p.meas_vt = (atoi(argv[25]) == 0 ? false : true); 
  
  // Lattice size 
  p.Nx = atoi(argv[32]);
  p.Ny = atoi(argv[33]);
  
  if(p.loop_max > std::min(p.Nx/2, p.Ny/2)) {
    cout << "Warning: requested Wilson loop max " << p.loop_max << " greater than ";
    cout << min(p.Nx/2, p.Ny/2) << ", truncating." << endl;
    p.loop_max = std::min(p.Nx/2, p.Ny/2);
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

  printParams(p);
  
  field<Complex> *gauge = new field<Complex>(p);
  gaussStart(gauge);  // hot start

  //Start simulation
  double time0 = -((double)clock());
  int iter_offset = 0;
  int iter = 0;
  cout << setprecision(16);
  
  leapfrogHMC *HMCStep = new leapfrogHMC(p);
  
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

    gettimeofday(&start, NULL);  
    for(iter=0; iter<2*p.therm; iter++){      
      //Perform HMC step
      accept = HMCStep->hmc(gauge, iter);
      gettimeofday(&total_end, NULL);  
      t_total = ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << t_total << " " << endl;                     //Time
    }
    gettimeofday(&end, NULL);  
    t_hmc += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
    iter_offset = 2*p.therm;
  }
    
  //Begin thermalised trajectories
  //---------------------------------------------------------------------
  for(iter=iter_offset; iter<p.iter_hmc + iter_offset; iter++){
    
    //Measure the topological charge at each step if trajectory is accepted
    //---------------------------------------------------------------------
    top = measTopCharge(gauge);
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
    
    //Perform Measurements
    //---------------------------------------------------------------------
    if((iter)%p.skip == 0 && iter > 2*p.therm) {
      
      count++; //Number of measurements taken
      
      //Plaquette action
      double plaq = measPlaq(gauge).real();
      plaqSum += plaq;

      //Dump simulation data to stdout
      gettimeofday(&total_end, NULL);  
      t_total = ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
      cout << fixed << setprecision(16) << iter << " ";   //Iteration
      cout << t_total << " ";                             //Time
      cout << plaqSum/count << " ";                       //Action
      cout << 1.0 - (double)top_stuck/(count*p.skip) << " " ;//P(Top transition)
      cout << HMCStep->exp_dH_ave/(count*p.skip) << " ";  //Average exp(-dH)
      cout << HMCStep->dH_ave/(count*p.skip) << " ";      //Average dH
      cout << (double)accepted/(count*p.skip) << " ";     //Acceptance
      cout << top_int << endl;                            //T charge
	
      //Dump simulation data to file
      name = "data/data/data"; //I cannot make bricks without clay!
      constructName(name, p);
      name += ".dat";	
      sprintf(fname, "%s", name.c_str());	
      fp = fopen(fname, "a");	
      fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e %.16e %d\n",
	      iter,
	      t_total,
	      plaqSum/count,
	      (double)top_stuck/(accepted),
	      HMCStep->exp_dH_ave/(count*p.skip),
	      HMCStep->dH_ave/(count*p.skip),
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
      
      //Creu<tz Ratios (for string tension)
      if(p.meas_wl) measWilsonLoops(gauge, plaq, iter);

      //Pion Correlation
      if(p.meas_pc) measPionCorrelation(gauge, iter);

      //Vacuum Trace
      //if(p.meas_vt) measVacuumTrace(gauge, top_old, iter, p);
      //-------------------------------------------------------------

      //Test deflation routines
      //-------------------------------------------------------------
      if(gauge->p.deflate) {
#if 1
	// Construct objects for an eigensolver
	//-------------------------------------
	eig_param_t eig_param;
	std::vector<field<Complex>*> kSpace;
	std::vector<Complex> evals;	
	prepareKrylovSpace(kSpace, evals, eig_param, gauge->p);
	
	// Compute a deflation space using IRAM
	iram(gauge, kSpace, evals, eig_param);
		
	// Test the block compression
	//---------------------------
	int Nx = gauge->p.Nx;
	int Ny = gauge->p.Ny;
	int Ns = 2;
	int blk_scheme[2] = {gauge->p.block_scheme[0], gauge->p.block_scheme[1]};
	int x_block_size = Nx/blk_scheme[0];
	int y_block_size = Ny/blk_scheme[1];
	int n_blocks = blk_scheme[0]*blk_scheme[1];
	int blk_size = Ns * x_block_size * y_block_size; // Complex elems per block
	int n_low = gauge->p.n_low;
	int n_conv = eig_param.n_conv;
	int n_deflate = eig_param.n_deflate;
	
	// Object to hold the block orthonormal low mode space
	std::vector<std::vector<Complex>> block_data_ortho(n_blocks, std::vector<Complex> (n_low * blk_size, 0.0));
	// Object to hold the projection coeffiecients of the high modes on the ow space
	std::vector<std::vector<Complex>> block_coef(n_blocks, std::vector<Complex> (n_low * n_conv, 0.0));

	gettimeofday(&start, NULL);

	// Krylov space
	std::vector<field<Complex>*> kSpace_recon(n_conv);
	for(int i=0; i<n_conv; i++) kSpace_recon[i] = new field<Complex>(gauge->p);
	// eigenvalues
	std::vector<Complex> evals_recon(n_conv);
	// Compress kSpace into block_data_ortho and block_coeffs...
	blockCompress(kSpace, block_data_ortho, block_coef, blk_scheme, n_low, n_conv);
	// ...then expand to into kSpace_recon test the quality
	blockExpand(kSpace_recon, block_data_ortho, block_coef, blk_scheme, n_low, n_conv);
	gettimeofday(&end, NULL);  
	
	// Compute the eigenvalues and residua using the reconstructed kSpace
	std::vector<double> resid(n_conv, 0.0);
	computeEvals(gauge, kSpace_recon, resid, evals_recon, n_conv);
	
	cout << "Compare eigenvalues and residua: " << endl;	
	for(int i=0; i<eig_param.n_conv; i++) printf("%d: %e %e \n", i, abs(evals[i].real() - evals_recon[i].real())/evals[i].real(), resid[i]);
	double delta_eval = 0.0;
	double delta_resid = 0.0;  
	for(int i=0; i<n_conv; i++) {
	  delta_eval += abs(evals[i].real() - evals_recon[i].real())/evals[i].real();
	  delta_resid += resid[i];
	}
	printf("<delta eval> = %e\n", delta_eval/n_conv);
	printf("<delta resid> = %e\n", delta_resid/n_conv);
	cout << endl;
	
	// Check compression ratio
	int pre = n_conv * 2 * Nx * Ny;
	int post= n_blocks * n_low * (blk_size + n_conv);
	cout << "Algorithmic compression: " << endl;
	cout << "Complex(double) elems pre = " << pre << " Complex(double) elems post = " << post << endl;
	cout << "Ratio: " << (100.0 * post)/pre << "% of original data " << endl;
	//cout << "Ratio2: " << 100*((1.0 * n_low)/n_conv + (1.0*n_low*n_blocks)/(Ns*Nx*Ny))<< "% of original data " << endl;
	cout << n_low << " low eigenvectors used " << endl;
	cout << (n_conv - n_low) << " high eigenvectors reconstructed " << endl;
	cout << "Compress/decompress time = " << ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 << endl;
	
	// Test deflation with true and reconstructed space
	//-------------------------------------------------
	field<Complex> *src = new field<Complex>(gauge->p);
	field<Complex> *sol = new field<Complex>(gauge->p);
	field<Complex> *check = new field<Complex>(gauge->p);

	// Populate src with rands
	gaussComplex(src);

	// Create inverter
	inverterCG *inv = new inverterCG(gauge->p);
	
	// Inversion with no deflation
	int undef_iter = inv->solve(sol, src, gauge);
	
	// Inversion with true deflation
	blas::zero(sol->data);
	int true_def_iter = inv->solve(sol, src, kSpace, evals, gauge);
	DdagDpsi(check, sol, gauge);
	blas::axpy(-1.0, src->data, check->data);
	
	cout << "Deflation efficacy: " << endl;
	cout << "Undeflated CG iter = " << undef_iter << endl;
	cout << "True Deflated CG iter   = " << true_def_iter << endl;
	cout << "Solution fidelity  = " << std::scientific << blas::norm2(check->data) << endl;

	// Inversion with true deflation
	blas::zero(sol->data);
	int recon_def_iter = inv->solve(sol, src, kSpace_recon, evals_recon, gauge);
	DdagDpsi(check, sol, gauge);
	blas::axpy(-1.0, src->data, check->data);
	
	cout << "Deflation efficacy: " << endl;
	cout << "Undeflated CG iter = " << undef_iter << endl;
	cout << "Recon Deflated CG iter = " << recon_def_iter << endl;
	cout << "Solution fidelity  = " << std::scientific << blas::norm2(check->data) << endl;
#endif
      }
      //-------------------------------------------------------------      
    }
    
    //Perform HMC step
    gettimeofday(&start, NULL);
    accept = HMCStep->hmc(gauge, iter);
    accepted += accept;
    gettimeofday(&end, NULL);  
    t_hmc += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
    //Checkpoint the gauge field?
    if((iter+1)%p.chkpt == 0) {	  
      name = "gauge/gauge";
      constructName(name, p);
      name += "_traj" + to_string(iter+1) + ".dat";
      writeGauge(gauge, name);
#ifdef HAVE_HDF5
      //hdf5Example();
#endif
    }
  }
  return 0;
}
//-------------------------------------------------------------------------------
