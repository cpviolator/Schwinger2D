#include "hmc.h"
#include "inverters.h"
#include "iram.h"

//2D HMC Routines
//---------------------------------------------------------------------

leapfrogHMC::leapfrogHMC(param_t param){
  inv = new inverterCG(param);

  guess_stack.reserve(4);
  for(int i=0; i<4; i++) {
    guess_stack.push_back(new field<Complex>(param));
  }
  phip = new field<Complex>(param);
  g3Dphi = new field<Complex>(param); 
  
};

bool leapfrogHMC::hmc_reversibility(field<Complex> *gauge, int iter) {

  int accept = 0;
  double H0, H1, H2;
  
  field<double> *mom = new field<double>(gauge->p);
  field<Complex> *phi = new field<Complex>(gauge->p);
  field<Complex> *chi = new field<Complex>(gauge->p);
  field<Complex> *gauge_old = new field<Complex>(gauge->p);  
  gauge_old->copy(gauge);
  
  // init mom[LX][LY][D]  <mom^2> = 1;
  gaussReal(mom); 
  
  if(gauge->p.dynamic) {    
    //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex(chi);
    //Create pseudo fermion field phi = D chi
    g3Dpsi(phi, chi, gauge);
  }

  // Forward trajectory
  H0 = measAction(mom, gauge, chi, false);
  trajectory(mom, gauge, phi, iter);
  H1 = measAction(mom, gauge, phi, true);

  // Reverse the trajectory
  gauge->p.tau *= -1.0;
  trajectory(mom, gauge, phi, iter);
  gauge->p.tau *= -1.0;
  H2 = measAction(mom, gauge, phi, true);

  // Metrics
  cout << "H0 = " << H0 << endl;
  cout << "H1 = " << H1 << endl;
  cout << "H2 = " << H2 << endl;
  cout << "H2 - H0 = " << std::scientific << H2 - H0 << endl;
  return (abs(H2 - H0) < 1e-4 ? true : false);
}

int leapfrogHMC::hmc(field<Complex> *gauge, int iter) {

  int accept = 0;
  double H = 0.0, H_old = 0.0;
  
  field<double> *mom = new field<double>(gauge->p);
  field<Complex> *phi = new field<Complex>(gauge->p);
  field<Complex> *chi = new field<Complex>(gauge->p);
  field<Complex> *gauge_old = new field<Complex>(gauge->p);  
  gauge_old->copy(gauge);
  
  // init mom[LX][LY][D]  <mom^2> = 1;
  gaussReal(mom); 
  
  if(gauge->p.dynamic) {    
    //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex(chi);
    //Create pseudo fermion field phi = D chi
    g3Dpsi(phi, chi, gauge);
  }

  // Reset the guess counter so that the previous guesses are discarded and new
  // ones constructed.
  guess_counter = 0;
  
  if (iter >= gauge->p.therm) H_old = measAction(mom, gauge, chi, false);
  trajectory(mom, gauge, phi, iter);
  if (iter >= gauge->p.therm) H = measAction(mom, gauge, phi, true);
  
  if (iter >= 2*gauge->p.therm) {      
    hmc_count++;
    exp_dH_ave += exp(-(H-H_old));
    dH_ave += (H-H_old);
  }

  // Metropolis accept/reject step
  if (iter >= gauge->p.therm) {    
    if ( drand48() > exp(-(H-H_old)) ) gauge->copy(gauge_old);
    else accept = 1;
  }
  
  return accept;
}

#if 0
void leapfrogHMC::trajectory(field<double> *mom, field<Complex> *gauge, field<Complex> *phi, int iter){
  
  double dtau = gauge->p.tau/gauge->p.n_step;
  double H = 0.0;
  //bool inspectrum_bool = true;
  bool inspectrum_bool = false;

  double ave_iter = 0;
  
  //gauge force (U field)
  field<double> *fU = new field<double>(gauge->p);
  //fermion force (D operator)
  field<double> *fD = new field<double>(gauge->p);

  // Compute the low spectrum
  //if(iter >= 2*gauge->p.therm && inspectrum_bool) inspectrum(gauge, iter);

  // Construct objects for an eigensolver
  //-------------------------------------
  eig_param_t eig_param;
  std::vector<field<Complex>*> kSpace;
  std::vector<field<Complex>*> kSpace_prior;
  std::vector<Complex> evals;
  std::vector<Complex> evals_prior;
  field<Complex> *gauge_prior;
  
  // Compute a deflation space using IRAM
  if(iter >= 2*gauge->p.therm && inspectrum_bool) {
    gauge_prior = new field<Complex>(gauge->p);
    gauge_prior->copy(gauge);
    prepareKrylovSpace(kSpace_prior, evals_prior, eig_param, gauge->p);
    //prepareKrylovSpace(kSpace, evals, eig_param, gauge->p);
    cout << "IRAM 0" << endl;
    iram(gauge_prior, kSpace_prior, evals_prior, eig_param);
  }

  //bool defl_update = false;
  
  // Start HMC trajectory
  //----------------------------------------------------------
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);
  update_mom(fU, fD, mom, 0.5*dtau);
  
  for(int k=1; k<gauge->p.n_step; k++) {

    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, dtau);
    
    if(iter >= 2*gauge->p.therm && inspectrum_bool) {
      cout << "IRAM " << k << endl;
      kspace_diff(gauge, kSpace, evals, kSpace_prior, evals_prior, eig_param, iter);
      gauge_prior->copy(gauge);
    }
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge);
    ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);    
    update_mom(fU, fD, mom, dtau);
  }

  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau);
  
  if(iter >= 2*gauge->p.therm && inspectrum_bool) {
    kspace_diff(gauge, kSpace, evals, kSpace_prior, evals_prior, eig_param, iter);
    gauge_prior->copy(gauge);
  }
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);
  update_mom(fU, fD, mom, 0.5*dtau);
  
  // HMC trajectory complete
  //----------------------------------------------------------
}
#endif

void leapfrogHMC::trajectory(field<double> *mom, field<Complex> *gauge, field<Complex> *phi, int iter){
  
  double dtau = gauge->p.tau/gauge->p.n_step;
  double H = 0.0;
  double ave_iter = 0;
  
  // gauge force (U field)
  field<double> *fU = new field<double>(gauge->p);
  // fermion force (D operator)
  field<double> *fD = new field<double>(gauge->p);
  // total force
  field<double> *f = new field<double>(gauge->p);

  // Construct objects for an eigensolver
  //-------------------------------------
  eig_param_t eig_param;
  std::vector<field<Complex>*> kSpace;
  std::vector<field<Complex>*> kSpace_prior;
  std::vector<Complex> evals;
  std::vector<Complex> evals_prior;
  field<Complex> *gauge_prior;

  double lambda = 1.0/6.0;
  double xi = 1.0/72.0;
  
  double lambda_dt = dtau*lambda;
  double dtauby2 = dtau / 2.0;
  double one_minus_2lambda_dt = (1-2*lambda)*dtau;
  double two_lambda_dt = lambda_dt*2;
  double xi_dtdt = 2*dtau*dtau*dtau*xi;  

  /*
  // Start LEAPFROG HMC trajectory
  //----------------------------------------------------------
  for(int k=1; k<=gauge->p.n_step; k++) {

    //Initial half step.
    //U_{k} = exp(i (dtau/2) P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, 0.5*dtau);
        
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge);
    ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);    
    update_mom(fU, fD, mom, dtau);

    //Final half step.
    //U_{k+1} = exp(i (dtau/2) P_{k+1/2}) * U_{k}
    update_gauge(gauge, mom, 0.5*dtau);
  }  
  // HMC trajectory complete
  //----------------------------------------------------------
  */
  
  // Start FGI trajectory
  //----------------------------------------------------------
  for(int k=1; k<=gauge->p.n_step; k++) {

    if(k == 1) {
      //forceU(fU, gauge);
      //update_mom(fU, mom, lambda_dt);    
      ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);
      update_mom(fD, mom, -lambda_dt);
    }

    innerFGI(fU, mom, gauge, dtauby2, gauge->p.inner_step);
    //update_gauge(gauge, mom, dtauby2);
    forceGradient(fU, fD, mom, phi, gauge, one_minus_2lambda_dt, xi_dtdt);
    innerFGI(fU, mom, gauge, dtauby2, gauge->p.inner_step);
    //update_gauge(gauge, mom, dtauby2);
    
    if(k == gauge->p.n_step) {
      //forceU(fU, gauge);
      //update_mom(fU, mom, lambda_dt);    
      ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);
      update_mom(fD, mom, -lambda_dt);
    } else {
      //forceU(fU, gauge);
      //update_mom(fU, mom, two_lambda_dt);    
      ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);
      update_mom(fD, mom, -two_lambda_dt);
    }
  }  
  
  // HMC trajectory complete
  //----------------------------------------------------------
}

void leapfrogHMC::forceGradient(field<double> *fU, field<double> *fD, field<double> *mom,
				field<Complex> *phi, field<Complex> *gauge, double one_minus_2lambda_dt, double xi_dtdt) {

  // Construct objects for an eigensolver
  //-------------------------------------
  std::vector<field<Complex>*> kSpace;
  std::vector<Complex> evals;
  
  // Copy the gauge and momentum, zero out the original momentum
  field<Complex> *gauge_copy = new field<Complex>(gauge->p);
  field<double> *mom_copy = new field<double>(mom->p);  
  gauge_copy->copy(gauge);
  mom_copy->copy(mom);
  blas::zero(mom->data);

  // Compute the forces
  //forceU(fU, gauge);
  //update_mom(fU, mom, xi_dtdt / one_minus_2lambda_dt);
  forceD(fD, phi, gauge, kSpace, evals, 0);
  update_mom(fD, mom, -xi_dtdt / one_minus_2lambda_dt);
  
  // Given the momentum kick (force), update the links to U'
  update_gauge(gauge, mom, 1.0);

  // Restore momentum
  mom->copy(mom_copy);

  // Add our kick to the momentum
  //forceU(fU, gauge);
  //update_mom(fU, mom, one_minus_2lambda_dt);
  forceD(fD, phi, gauge, kSpace, evals, 0);
  update_mom(fD, mom, -one_minus_2lambda_dt);
  
  // Restore the gauge field
  gauge->copy(gauge_copy);

  delete gauge_copy;
  delete mom_copy;
}

void leapfrogHMC::forceGradient(field<double> *fU, field<double> *mom, field<Complex> *gauge, double one_minus_2lambda_dt, double xi_dtdt) {

  // Construct objects for an eigensolver
  //-------------------------------------
  std::vector<field<Complex>*> kSpace;
  std::vector<Complex> evals;
  
  // Copy the gauge and momentum, zero out the original momentum
  field<Complex> *gauge_copy = new field<Complex>(gauge->p);
  field<double> *mom_copy = new field<double>(mom->p);  
  gauge_copy->copy(gauge);
  mom_copy->copy(mom);
  blas::zero(mom->data);

  // Compute the forces
  forceU(fU, gauge);
  update_mom(fU, mom, xi_dtdt / one_minus_2lambda_dt);
  
  // Given the momentum kick (force), update the links to U'
  update_gauge(gauge, mom, 1.0);

  // Restore momentum
  mom->copy(mom_copy);

  // Add our kick to the momentum
  forceU(fU, gauge);
  update_mom(fU, mom, one_minus_2lambda_dt);
  
  // Restore the gauge field
  gauge->copy(gauge_copy);

  delete gauge_copy;
  delete mom_copy;
}


void leapfrogHMC::innerFGI(field<double> *fU, field<double> *mom, field<Complex> *gauge, double tau, int steps) {

  double lambda = 1.0/6.0;
  double xi = 1.0/72.0;

  double dtau = tau / steps;
  double lambda_dt = dtau*lambda;
  double dtauby2 = dtau / 2.0;
  double one_minus_2lambda_dt = (1-2*lambda)*dtau;
  double two_lambda_dt = lambda_dt*2;
  double xi_dtdt = 2*dtau*dtau*dtau*xi;

  for(int k=1; k <= steps; k++){

    if(k == 1) {
      forceU(fU, gauge);
      update_mom(fU, mom, lambda_dt);
    }

    update_gauge(gauge, mom, dtauby2);
    forceGradient(fU, mom, gauge, one_minus_2lambda_dt, xi_dtdt);
    update_gauge(gauge, mom, dtauby2);

    forceU(fU, gauge);
    if(k == steps) update_mom(fU, mom, lambda_dt);
    else           update_mom(fU, mom, two_lambda_dt);      
    
  }
}

void leapfrogHMC::forceU(field<double> *fU, field<Complex> *gauge) {
  
  Complex plaq0 = 0.0, plaq = 0.0;
  double temp = 0.0;
  int xp1, xm1, yp1, ym1;

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  double beta = gauge->p.beta;
  for(int x=0; x<Nx; x++) {
    xp1 = (x+1)%Nx;
    xm1 = (x-1+Nx)%Nx;
    
    for(int y=0; y<Ny; y++) {

      yp1 = (y+1)%Ny;
      ym1 = (y-1+Ny)%Ny;
      
      plaq0 = gauge->read(x,y,0)*gauge->read(xp1,y,1)*conj(gauge->read(x,yp1,0))*conj(gauge->read(x,y,1));
      plaq =  gauge->read(x,ym1,0)*gauge->read(xp1,ym1,1)*conj(gauge->read(x,y,0))*conj(gauge->read(x,ym1,1));
      temp = beta*(imag(plaq0) - imag(plaq));
      fU->write(x,y,0, temp);
      
      plaq = gauge->read(x,y,1)*conj(gauge->read(xm1,yp1,0))*conj(gauge->read(xm1,y,1))*gauge->read(xm1,y,0);
      temp = beta*(imag(plaq) - imag(plaq0));
      fU->write(x,y,1, temp);
    }
  }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void leapfrogHMC::update_mom(field<double> *fU, field<double> *fD, field<double> *mom, double dtau){
  
  int Nx = fU->p.Nx;
  int Ny = fU->p.Ny;
  double temp = 0.0;
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
	temp = mom->read(x,y,mu) - (fU->read(x,y,mu) - fD->read(x,y,mu))*dtau;
	mom->write(x,y,mu, temp);
      }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (f)
void leapfrogHMC::update_mom(field<double> *f, field<double> *mom, double dtau){
  
  int Nx = f->p.Nx;
  int Ny = f->p.Ny;
  double temp = 0.0;
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
	temp = mom->read(x,y,mu) - dtau*(f->read(x,y,mu));
	mom->write(x,y,mu, temp);
      }
}


//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
// MILC: update_u
void leapfrogHMC::update_gauge(field<Complex> *gauge, field<double> *mom, double dtau){
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  Complex temp = 0.0;  
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
	temp = gauge->read(x,y,mu) * polar(1.0, mom->read(x,y,mu) * dtau);
	gauge->write(x,y,mu, temp);
      }
}

void leapfrogHMC::update_deflation(field<Complex> *gauge, field<Complex> *gauge_prior, std::vector<field<Complex>*> &kSpace, std::vector<Complex> &evals){
  
  int n_vecs = evals.size();
  std::vector<field<Complex>*> kSpace_copy;
  kSpace_copy.reserve(n_vecs);
  for(int i=0; i<n_vecs; i++) {
    kSpace_copy.push_back(new field<Complex>(gauge->p));
    kSpace_copy[i]->copy(kSpace[i]);
    //cout << " Norm of kSpace[" << i << "] = " << blas::norm(kSpace_copy[i]->data) << endl;
  }
  
  // Compute rotation
  Eigen::MatrixXcd epsilon = MatrixXcd::Zero(n_vecs, n_vecs);
  field<Complex> *temp1 = new field<Complex>(gauge->p);
  field<Complex> *temp2 = new field<Complex>(gauge->p);
  for(int i=0; i<n_vecs; i++) {
    DdagDpsi(temp1, kSpace[i], gauge);
    DdagDpsi(temp2, kSpace[i], gauge_prior);
    blas::axpy(-1.0, temp2->data, temp1->data);
    for(int j=0; j<n_vecs; j++) {
      epsilon(i,j) = blas::cDotProd(kSpace[j]->data, temp1->data);
      //cout << "epsilon("<<i<<","<<j<<") = " << epsilon(i,j) << endl;
    }
  }

  // Rotate vectors
  for(int i=0; i<n_vecs; i++) {
    for(int j=0; j<n_vecs && j != i; j++) {
      blas::caxpy(epsilon(i,j) / (evals[i].real() - evals[j].real()), kSpace[j]->data, kSpace_copy[i]->data);
    }
  }

  // Peturb eigenvalues, restore krylov space
  for(int i=0; i<n_vecs; i++) {
    cout << evals[i] << " + " << epsilon(i,i) << endl;
    evals[i] += epsilon(i,i);
    kSpace[i]->copy(kSpace_copy[i]);
    blas::ax(1.0/blas::norm(kSpace[i]->data), kSpace[i]->data);
  }

  // Test fidelity
  for (int i = 0; i < n_vecs; i++) {
    
    // r = A * v_i
    DdagDpsi(temp1, kSpace[i], gauge);
    
    // lambda_i = v_i^dag A v_i / (v_i^dag * v_i)
    // Measure ||lambda_i*v_i - A*v_i||
    Complex n_unit(-1.0, 0.0);
    blas::caxpby(blas::cDotProd(kSpace[i]->data, temp1->data), kSpace[i]->data, n_unit, temp1->data);
    cout << "residua = " << blas::norm(temp1->data) << endl;
  }
}

void leapfrogHMC::kspace_diff(field<Complex> *gauge, std::vector<field<Complex>*> &kSpace, std::vector<Complex> &evals,
			      std::vector<field<Complex>*> &kSpace_prior, std::vector<Complex> &evals_prior, eig_param_t &param, int iter){
  
  int n_vecs = evals.size();

  // Compute the eigendecomposition of the new gauge field
  prepareKrylovSpace(kSpace, evals, eig_param, gauge->p);
  iram(gauge, kSpace, evals, eig_param);
  writeEigenData(kSpace, evals, eig_param, iter);
  
  // Compute orthonormality matrix
  for(int i=0; i<n_vecs; i++) {
    for(int j=0; j<n_vecs; j++) {
      cout << abs(blas::cDotProd(kSpace[i]->data, kSpace_prior[j]->data)) << " ";
    }
    cout << endl;
  }
  
  // Copy the new space to the prior space
  for(int i=0; i<n_vecs; i++) {      
    evals_prior[i] = evals[i];
    kSpace_prior[i]->copy(kSpace[i]);  
  }
}


// Optimise this to operate only on a single parity of sites.
int leapfrogHMC::forceD(field<double> *fD, field<Complex> *phi, field<Complex> *gauge,
			std::vector<field<Complex>*> &kSpace, std::vector<Complex> &evals, int iter)
{
  int cg_iter = 0;
  if(gauge->p.dynamic == true) {

    blas::zero(fD->data);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1)
    //  phi = (D^-1 * g3 * D^-1 g3) phi.

    // Inspect guess
    //for(int i=0; i<10; i++) phip->print(i);
    if(iter < 2*gauge->p.therm || !gauge->p.deflate) {
      cg_iter += inv->solve(phip, phi, gauge);
    } else {
      cg_iter += inv->solve(phip, phi, kSpace, evals, gauge);
    }
    
    //cg_iter += inv->solve(phip, phi, gauge);
    //for(int i=0; i<10; i++) phip->print(i);
    
    //g3Dphi = g3D * phip
    g3Dpsi(g3Dphi, phip, gauge);

    Complex gauge_px, gauge_py, phip0, phip0_px, phip0_py,
      phip1, phip1_px, phip1_py,
      g3Dphi0, g3Dphi0_px, g3Dphi0_py,
      g3Dphi1, g3Dphi1_px, g3Dphi1_py;
    
    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
#pragma omp parallel for
    for(int x=0; x<Nx; x++) {
      int xp1 = (x+1)%Nx;
      for(int y=0; y<Ny; y++) {
	int yp1 = (y+1)%Ny;

	double temp = 0.0;
	
	/*
	// Collect all data
	gauge_px = gauge->read(x,y,0);
	gauge_py = gauge->read(x,y,1);
	
	phip0 =    conj(phip->read(x,  y,0));
	phip0_px = conj(phip->read(xp1,y,0));
	phip0_py = conj(phip->read(x,yp1,0));
	
	phip1 =    conj(phip->read(x,  y,1));
	phip1_px = conj(phip->read(xp1,y,1));
	phip1_py = conj(phip->read(x,yp1,1));
	
	g3Dphi0 =    g3Dphi->read(x,  y,0);
	g3Dphi0_px = g3Dphi->read(xp1,y,0);
	g3Dphi0_py = g3Dphi->read(x,yp1,0);
	
	g3Dphi1 =    g3Dphi->read(x,  y,1);
	g3Dphi1_px = g3Dphi->read(xp1,y,1);
	g3Dphi1_py = g3Dphi->read(x,yp1,1);

	
	// Apply antiperiodic boundary conditions in y direction
	if (yp1 == 0) {
	  //phip0_py *= -1.0;
	  //phip1_py *= -1.0;
	  //g3Dphi0_py *= -1.0;
	  //g3Dphi1_py *= -1.0;
	}

	//mu = 0
	//upper
	// | r  1 |
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |
	temp = real(I*((conj(gauge_px) *
			(phip0_px * (r*g3Dphi0 +   g3Dphi1) -
			 phip1_px * (  g3Dphi0 + r*g3Dphi1)))
		       - (gauge_px *
			  (phip0 * (r*g3Dphi0_px -   g3Dphi1_px) +
			   phip1 * (  g3Dphi0_px - r*g3Dphi1_px)))));

	fD->write(x,y,0,temp);
	
	//mu = 1
	//upper
	// | r -i |
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	temp = real(I*((conj(gauge_py) *
			(phip0_py * (r*g3Dphi0 - I*g3Dphi1) -
			 phip1_py * (I*g3Dphi0 + r*g3Dphi1)))
		       - (gauge_py *
			  (phip0 * (r*g3Dphi0_py + I*g3Dphi1_py) +
			   phip1 * (I*g3Dphi0_py - r*g3Dphi1_py)))));
	
	fD->write(x,y,1,temp);
	*/
	
	//mu = 0
	//upper
	// | r  1 | 
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |	
	temp = real(I*((conj(gauge->read(x,y,0)) *
			(conj(phip->read(xp1,y,0)) * (r*g3Dphi->read(x,y,0) +   g3Dphi->read(x,y,1)) -
			 conj(phip->read(xp1,y,1)) * (  g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
		       -
		       (gauge->read(x,y,0) *
			(conj(phip->read(x,y,0)) * (r*g3Dphi->read(xp1,y,0) -   g3Dphi->read(xp1,y,1)) +
			 conj(phip->read(x,y,1)) * (  g3Dphi->read(xp1,y,0) - r*g3Dphi->read(xp1,y,1))))
		       )
		    );
	
	fD->write(x,y,0,temp);
	
	//mu = 1
	//upper
	// | r -i | 
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	temp = real(I*((conj(gauge->read(x,y,1)) *
			(conj(phip->read(x,yp1,0)) * (r*g3Dphi->read(x,y,0) - I*g3Dphi->read(x,y,1)) -
			 conj(phip->read(x,yp1,1)) * (I*g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
		       -			       
		       (gauge->read(x,y,1) *
			(conj(phip->read(x,y,0)) * (r*g3Dphi->read(x,yp1,0) + I*g3Dphi->read(x,yp1,1)) +
			 conj(phip->read(x,y,1)) * (I*g3Dphi->read(x,yp1,0) - r*g3Dphi->read(x,yp1,1))))
		       )
		    );
	
	fD->write(x,y,1,temp);
	
      }
    }
    /*
    // phip and g3phi have now been read and may be used for the next step
    guess_stack[guess_counter]->copy(phip);
    if(guess_counter == 1) {
      //cout << "Shifting " << endl;
      // guess_delta = guess_{n} - guess_{n-1}
      blas::axpy(-1.0, guess_stack[0]->data, guess_stack[1]->data, guess_stack[2]->data);

      // guess_{n+1} = solution_{n} + guess_delta
      blas::axpy(1.0, guess_stack[2]->data, phip->data);

      // shift the stack,
      guess_stack[0]->copy(guess_stack[1]);
      
    } else {
      //cout << "Incrementing " << endl;
      guess_counter++;
    } 
    */
    //blas::zero(phip->data);
  }
  return cg_iter;
}



//----------------------------------------------------------------------------------

leapfrogHMC::~leapfrogHMC() {
  
  evals0.resize(0);
  for (int i=0; i<kSpace0.size(); i++) {
    delete kSpace0[i];
  }
  
  evals1.resize(0);
  for (int i=0; i<kSpace1.size(); i++) {
    delete kSpace1[i];
  }
  
  evals_delta.resize(0);
  for (int i=0; i<kSpace_delta.size(); i++) {
    delete kSpace_delta[i];
  }
  
  evals_prediction.resize(0);
  for (int i=0; i<kSpace_prediction.size(); i++) {
    delete kSpace_prediction[i];
  }
  
};
