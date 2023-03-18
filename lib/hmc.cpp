#include "hmc.h"
#include "inverters.h"
#include "io.h"

//2D HMC Routines
//---------------------------------------------------------------------
HMC::HMC(Param param) {

  guess_stack.reserve(10);  
  for(int i=0; i<10; i++) guess_stack.push_back(new field<Complex>(param));

  phip = new field<Complex>(param);
  g3Dphi = new field<Complex>(param);

  if(param.flavours == 3) {
#ifdef ENABLE_ALG_REMEZ
    int n = param.pfe_degree; // The degree of the numerator polynomial
    int d = param.pfe_degree; // The degree of the denominator polynomial
    int y = 1;  // The numerator of the exponent
    int z = 4;  // The denominator of the exponent
    int precision = param.pfe_prec; // The precision that gmp uses
    double lambda_low = 0.0001, lambda_high = 32; // The bounds of the approximation
    
    // Instantiate the Remez class
    remez = new AlgRemez(lambda_low, lambda_high, precision);
    
    // The partial fraction expansion takes the form 
    // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
    
    // Heatbath PFE
    heatbath_pfe.inv_res.resize(n);
    heatbath_pfe.inv_pole.resize(d);
    heatbath_pfe.res.resize(n);
    heatbath_pfe.pole.resize(d);

    string heatbath_pfe_name;    
    heatbath_pfe_name += "remez_pfe_wisdom_n" + to_string(n) + "_d" + to_string(d) + "_y" + to_string(y) + "_z"+ to_string(z) + "_prec" + to_string(precision) + "_Llow" + to_string(lambda_low) + "_Lhigh" + to_string(lambda_high) + ".dat";
    
    // Read or generate the required approximation
    if(!readPFE(heatbath_pfe, heatbath_pfe_name)) {
      remez->generateApprox(n,d,y,z);
      // Find the partial fraction expansion of the approximation 
      // to the function x^{y/z} (this only works currently for 
      // the special case that n = d)
      remez->getIPFE(heatbath_pfe.inv_res.data(), heatbath_pfe.inv_pole.data(), &heatbath_pfe.inv_norm);
      remez->getPFE(heatbath_pfe.res.data(), heatbath_pfe.pole.data(), &heatbath_pfe.norm);
      writePFE(heatbath_pfe, heatbath_pfe_name);
    }
    printf("IPFE:\nalpha[0] = %18.16e\n", heatbath_pfe.inv_norm);
    for (int i = 0; i < n; i++) {
      printf("alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	     i+1, heatbath_pfe.inv_res[i], i+1, heatbath_pfe.inv_pole[i]);
    }
    printf("PFE:\nalpha[0] = %18.16e\n", heatbath_pfe.norm);
    for (int i = 0; i < n; i++) {
      printf("alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	     i+1, heatbath_pfe.res[i], i+1, heatbath_pfe.pole[i]);
    }
    
    // Force PFE
    force_pfe.inv_res.resize(n);
    force_pfe.inv_pole.resize(d);
    force_pfe.res.resize(n);
    force_pfe.pole.resize(d);

    y = 1;
    z = 2;
    // Read or generate the required approximation
    string force_pfe_name;    
    force_pfe_name += "remez_pfe_wisdom_n" + to_string(n) + "_d" + to_string(d) + "_y" + to_string(y) + "_z"+ to_string(z) + "_prec" + to_string(precision) + "_Llow" + to_string(lambda_low) + "_Lhigh" + to_string(lambda_high) + ".dat";

    if(!readPFE(force_pfe, force_pfe_name)) {
      remez->generateApprox(n,d,y,z);
      // Find the partial fraction expansion of the approximation 
      // to the function x^{y/z} (this only works currently for 
      // the special case that n = d)
      remez->getIPFE(force_pfe.inv_res.data(), force_pfe.inv_pole.data(), &force_pfe.inv_norm);
      remez->getPFE(force_pfe.res.data(), force_pfe.pole.data(), &force_pfe.norm);
      writePFE(force_pfe, force_pfe_name);
    }
    printf("IPFE:\nalpha[0] = %18.16e\n", force_pfe.inv_norm);
    for (int i = 0; i < n; i++) {
      printf("alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	     i+1, force_pfe.inv_res[i], i+1, force_pfe.inv_pole[i]);
    }
    printf("PFE:\nalpha[0] = %18.16e\n", force_pfe.norm);
    for (int i = 0; i < n; i++) {
      printf("alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	     i+1, force_pfe.res[i], i+1, force_pfe.pole[i]);
    }
#else
    cout << "Error: AlgRemez not installed. Please recompile with ENABLE_ALG_REMEZ = ON: flavours passed = " << param.flavours << endl;
    exit(0);
#endif
  } else if(!(param.flavours == 2 || param.flavours == 0)) {
    cout << "Error: flavours must be 0, 2, or 3: flavours passed = " << param.flavours << endl;
    exit(0);    
  }
};

bool HMC::hmc_reversibility(field<Complex> *gauge, int iter) {

  int accept = 0;
  double H0, H1, H2;

  field<double> *mom = new field<double>(gauge->p);
  std::vector<field<Complex>*> phi;
  std::vector<field<Complex>*> chi;
  phi.reserve(2); chi.reserve(2);
  for(int i=0; i<2; i++) {
    phi.push_back(new field<Complex>(gauge->p));
    chi.push_back(new field<Complex>(gauge->p));
  }
  field<Complex> *gauge_old = new field<Complex>(gauge->p);  
  blas::copy(gauge_old, gauge);
  
  // init mom[LX][LY][D]  <mom^2> = 1;
  gaussReal(mom); 

  // Create CG solver
  inv = new inverterCG(gauge->p);
  
  if(gauge->p.flavours > 0) {    
    //Create gaussian distributed fermion field chi
    gaussComplex(chi[0]);
    if(gauge->p.flavours == 3) gaussComplex(chi[1]);
    
    //Create pseudo fermion field phi = g3Ddag chi
    g3Dpsi(phi[0], chi[0], gauge);
    if(gauge->p.flavours == 3) {
      //Create a rational pseudo fermion field phi = (g3 D g3 D)^-1/4 chi
      std::vector<field<Complex>*> phi_arr;
      for(int i=0; i<gauge->p.pfe_degree; i++) phi_arr.push_back(new field<Complex>(gauge->p));    
      inv->solveMulti(phi_arr, chi[1], gauge, heatbath_pfe.pole);
      blas::zero(phi[1]);
      // Accumulate the norm
      blas::axpy(heatbath_pfe.norm, chi[1], phi[1]);
      // Continue accumulation
      for(int i=0; i<gauge->p.pfe_degree; i++) blas::axpy(heatbath_pfe.res[i], phi_arr[i], phi[1]);
      for(int i=0; i<gauge->p.pfe_degree; i++) delete phi_arr[i];        
    }
  }

  // Forward trajectory
  // Temporarily swich off deflation and inspection
  inv->switchOffEigensolver();  
  H0 = measAction(mom, gauge, phi, heatbath_pfe);
  inv->switchOnEigensolver();
  
  trajectory(mom, gauge, phi, iter);

  inv->switchOffEigensolver();  
  H1 = measAction(mom, gauge, phi, heatbath_pfe);  
  inv->switchOnEigensolver();
  
  // Reverse the trajectory
  gauge->p.tau *= -1.0;
  trajectory(mom, gauge, phi, iter);
  gauge->p.tau *= -1.0;
  H2 = measAction(mom, gauge, phi, heatbath_pfe);

  // Metrics
  cout << "H0 = " << H0 << endl;
  cout << "H1 = " << H1 << endl;
  cout << "H2 = " << H2 << endl;
  cout << "H2 - H0 = " << std::scientific << H2 - H0 << endl;
  printf("(H2 - H0)/H0 percent error = %.12e\n", 100*abs(H2 - H0) / H0);

  delete mom;
  delete phi[0]; delete phi[1];
  delete chi[0]; delete chi[1];
  delete gauge_old;
  return (abs(H2 - H0) < 1e-4 ? true : false);
}

int HMC::hmc(field<Complex> *gauge, int iter) {

  int accept = 0;
  double H = 0.0, H_old = 0.0;
  
  field<double> *mom = new field<double>(gauge->p);
  std::vector<field<Complex>*> phi;
  std::vector<field<Complex>*> chi;
  phi.reserve(2); chi.reserve(2);
  for(int i=0; i<2; i++) {
    phi.push_back(new field<Complex>(gauge->p));
    chi.push_back(new field<Complex>(gauge->p));
  }
  
  field<Complex> *gauge_old = new field<Complex>(gauge->p);  
  blas::copy(gauge_old, gauge);
  
  // init momentum
  gaussReal(mom);

  // Create CG solver
  inv = new inverterCG(gauge->p);
  
  if(gauge->p.flavours > 0) {    
    //Create gaussian distributed fermion field chi
    gaussComplex(chi[0]);
    if(gauge->p.flavours == 3) gaussComplex(chi[1]);
    
    //Create pseudo fermion field phi = g3Ddag chi
    g3Dpsi(phi[0], chi[0], gauge);
    if(gauge->p.flavours == 3) {
      //Create a rational pseudo fermion field phi = (g3 D g3 D)^-1/4 chi
      std::vector<field<Complex>*> phi_arr;
      for(int i=0; i<gauge->p.pfe_degree; i++) phi_arr.push_back(new field<Complex>(gauge->p));    
      inv->solveMulti(phi_arr, chi[1], gauge, heatbath_pfe.pole);
      blas::zero(phi[1]);
      // Accumulate the norm
      blas::axpy(heatbath_pfe.norm, chi[1], phi[1]);
      // Continue accumulation
      for(int i=0; i<gauge->p.pfe_degree; i++) blas::axpy(heatbath_pfe.res[i], phi_arr[i], phi[1]);
      for(int i=0; i<gauge->p.pfe_degree; i++) delete phi_arr[i];        
    }
  }

  if(iter >= gauge->p.therm) {    
    // Temporarily swich off deflation and inspection
    inv->switchOffEigensolver();      
    // H_old = P^2 + S(U) + <chi|chi>
    H_old = measAction(mom, gauge, phi, heatbath_pfe);
    inv->switchOnEigensolver();
  }
  
  // Perfrom trajectory
  trajectory(mom, gauge, phi, iter);
  if ( !(gauge->p.sampler == S_HMC) && (gauge->p.sampler == S_MCHMC)) langevin_noise(mom,gauge);
      
  
  if(iter >= gauge->p.therm) {
    // H_evolved = P^2 + S(U) + <phi| (Ddag D)^-1 |phi>
    inv->switchOffEigensolver();  
    H = measAction(mom, gauge, phi, heatbath_pfe);
    inv->switchOnEigensolver();  
  }
  
  if (iter >= gauge->p.therm) {      
    hmc_count++;
  }
  
  if (gauge->p.sampler == S_HMC ) exp_dH = exp(-(H-H_old));
  dH = (H-H_old);
    std::cout<<dH<<"\n";
  // Metropolis accept/reject step
  if (iter >= gauge->p.therm) {    
    if ( ( drand48() > exp(-(H-H_old)) ) && !(gauge->p.sampler == S_MCHMC ) ) {

      // Reject the configuration
      blas::copy(gauge, gauge_old);
    }
    else {
      // Accept this configuration
      accept = 1;
      // If a trajectory is rejected, the action difference
      // is zero, so we only add when the gauge is accepted.
      exp_dH_ave += exp_dH;
      dH_ave += dH;
    }
  }

  delete mom;
  delete phi[0]; delete phi[1];
  delete chi[0]; delete chi[1];
  delete gauge_old;
  delete inv;
  
  return accept;
}

field<double>* HMC::computeFermionForce(field<Complex> *gauge, std::vector<field<Complex>*> &phi) {
  
  field<double> *fD = new field<double>(gauge->p);
  field<double> *fD_rat;
  
  forceD(fD, phi[0], gauge);
  if(gauge->p.flavours == 3) {
    fD_rat = new field<double>(gauge->p);
    forceMultiD(fD_rat, phi[1], gauge);
    blas::axpy(1.0, fD_rat, fD);
  }
  return fD;
}

void HMC::trajectory(field<double> *mom, field<Complex> *gauge, std::vector<field<Complex>*> &phi, int iter){
  
  double dtau = gauge->p.tau/gauge->p.n_step;
  double H = 0.0;
  double ave_iter = 0;
  
  // gauge force (U field)
  field<double> *fU;
  
  // fermion force (D operator)
  field<double> *fD;

  double lambda = 1.0/6.0;
  double xi = 1.0/72.0;  
  double lambda_dt = dtau*lambda;
  double dtauby2 = dtau / 2.0;
  double one_minus_2lambda_dt = (1-2*lambda)*dtau;
  double two_lambda_dt = lambda_dt*2;
  double xi_dtdt = 2*dtau*dtau*dtau*xi;
  double inner_step = gauge->p.inner_step;
  
  switch(gauge->p.integrator) {
  case LEAPFROG:
    // Start LEAPFROG HMC trajectory
    //----------------------------------------------------------
    for(int k=0; k<gauge->p.n_step; k++) {
      
      //Initial half step.
      //U_{k} = exp(i (dtau/2) P_{k-1/2}) * U_{k-1}
      update_gauge(gauge, mom, dtauby2);
      
      //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
      fU = computeGaugeForce(gauge);
      fD = computeFermionForce(gauge, phi);
      blas::axpy(-1.0, fD, fU);
      if (gauge->p.sampler == S_HMC) 
          update_mom_s(fU, mom, dtau,gauge->p.sampler);
      else if (gauge->p.sampler == S_MCHMC) 
          update_mom_mchmc(fU, mom, dtau);
      
      //Final half step.
      //U_{k+1} = exp(i (dtau/2) P_{k+1/2}) * U_{k}
      update_gauge(gauge, mom, dtauby2);
    }    
    break;

  case FGI:  
    // Start FGI trajectory
    //----------------------------------------------------------
    for(int k=0; k<gauge->p.n_step; k++) {

#if 0
      // QPQPQ: less efficient
      innerFGI(mom, gauge, dtauby2/2, inner_step);
      forceGradient(mom, phi, gauge, one_minus_2lambda_dt/2, xi_dtdt/2);
      innerFGI(mom, gauge, dtauby2/2, inner_step);

      fD = computeFermionForce(gauge, phi);
      update_mom_s(fD, mom, -two_lambda_dt,gauge->p.sampler);
      
      innerFGI(mom, gauge, dtauby2/2, inner_step);
      forceGradient(mom, phi, gauge, one_minus_2lambda_dt/2, xi_dtdt/2);
      innerFGI(mom, gauge, dtauby2/2, inner_step);
#else
      //PQPQP: more efficient
      if(k == 0) {
	fD = computeFermionForce(gauge, phi);
	update_mom_s(fD, mom, -lambda_dt, gauge->p.sampler);
      }
      
      innerFGI(mom, gauge, dtauby2, inner_step);
      forceGradient(mom, phi, gauge, one_minus_2lambda_dt, xi_dtdt);
      innerFGI(mom, gauge, dtauby2, inner_step);
      
      fD = computeFermionForce(gauge, phi);
      update_mom_s(fD, mom, k == gauge->p.n_step - 1 ? -lambda_dt : -two_lambda_dt,gauge->p.sampler);
#endif
    }    
    break;
  default: cout << "Error: unknown integrator type" << endl; exit(0);
  }
  // HMC trajectory complete
  //----------------------------------------------------------
}

void HMC::forceGradient(field<double> *mom, std::vector<field<Complex>*> phi, field<Complex> *gauge,
			double one_minus_2lambda_dt, double xi_dtdt) {
  
  // Copy the gauge and momentum, zero out the original momentum
  field<Complex> *gauge_copy = new field<Complex>(gauge->p);
  field<double> *mom_copy = new field<double>(mom->p);
  field<double> *fD;
  
  blas::copy(gauge_copy, gauge);
  blas::copy(mom_copy, mom);
  blas::zero(mom);

  fD = computeFermionForce(gauge, phi);
  update_mom_s(fD, mom, -xi_dtdt / one_minus_2lambda_dt, gauge->p.sampler);
  
  // Given the momentum kick (force), update the links to U'
  update_gauge(gauge, mom, 1.0);

  // Restore momentum
  blas::copy(mom, mom_copy);
  
  // Add our kick to the momentum
  fD = computeFermionForce(gauge, phi);
  update_mom_s(fD, mom, -one_minus_2lambda_dt, gauge->p.sampler);
  
  // Restore the gauge field
  blas::copy(gauge, gauge_copy);
  
  delete gauge_copy;
  delete mom_copy;
}

void HMC::forceGradient(field<double> *mom, field<Complex> *gauge,
			double one_minus_2lambda_dt, double xi_dtdt) {

  // Copy the gauge and momentum, zero out the original momentum
  field<Complex> *gauge_copy = new field<Complex>(gauge->p);
  field<double> *mom_copy = new field<double>(mom->p);
  field<double> *fU;
  blas::copy(gauge_copy, gauge);
  blas::copy(mom_copy, mom);
  blas::zero(mom);

  // Compute the forces
  fU = computeGaugeForce(gauge);
  update_mom_s(fU, mom, xi_dtdt / one_minus_2lambda_dt, gauge->p.sampler);
  
  // Given the momentum kick (force), update the links to U'
  update_gauge(gauge, mom, 1.0);

  // Restore momentum
  blas::copy(mom, mom_copy);

  // Add our kick to the momentum
  fU = computeGaugeForce(gauge);
  update_mom_s(fU, mom, one_minus_2lambda_dt, gauge->p.sampler);
  
  // Restore the gauge field
  blas::copy(gauge, gauge_copy);
  
  delete gauge_copy;
  delete mom_copy;
}


void HMC::innerFGI(field<double> *mom, field<Complex> *gauge, double tau, int steps) {

  double lambda = 1.0/6.0;
  double xi = 1.0/72.0;

  double dtau = tau / steps;
  double lambda_dt = dtau*lambda;
  double dtauby2 = dtau / 2.0;
  double one_minus_2lambda_dt = (1-2*lambda)*dtau;
  double two_lambda_dt = lambda_dt*2;
  double xi_dtdt = 2*dtau*dtau*dtau*xi;

  field<double> *fU;
  
  for(int k=0; k < steps; k++){

    if(k == 0) {
      fU = computeGaugeForce(gauge);
      update_mom_s(fU, mom, lambda_dt, gauge->p.sampler);
    }

    update_gauge(gauge, mom, dtauby2);
    forceGradient(mom, gauge, one_minus_2lambda_dt, xi_dtdt);
    update_gauge(gauge, mom, dtauby2);

    fU = computeGaugeForce(gauge);
    update_mom_s(fU, mom, k == steps - 1 ? lambda_dt : two_lambda_dt, gauge->p.sampler);
    
  }
}

field<double>* HMC::computeGaugeForce(field<Complex> *gauge) {

  field<double> *fU = new field<double>(gauge->p);
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
  return fU;
}
//P_{k+1/2} = P_{k-1/2} - dtau * (f)
void HMC::update_mom(field<double> *f, field<double> *mom, double dtau){
  
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



//P_{k+1/2} = P_{k-1/2} - dtau * (f)
void HMC::update_mom_s(field<double> *f, field<double> *mom, double dtau, Sampler sampler){
  
  int Nx = f->p.Nx;
  int Ny = f->p.Ny;
  int d = Nx * Ny * 2;
  double temp = 0.0;
   
  if (sampler == S_HMC) {  
      for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
          for(int mu=0; mu<2; mu++) {
        temp = mom->read(x,y,mu) - dtau*(f->read(x,y,mu));
        mom->write(x,y,mu, temp);
          }
  }
  else if (sampler == S_MCHMC) {
      
      int f_norm = sqrt(blas::norm2(f)); 
      double ue = -1.0 * blas::DotProd(f,mom)/f_norm;
      double sh = sinh(dtau * f_norm / (d-1));
      double ch = cosh(dtau * f_norm / (d-1));
      double th = tanh(dtau * f_norm / (d-1));
      double delta_r = log(ch) + log1p(ue*th);
      
//                   double g_norm = sqrt(blas::norm2(f)); 

//             //update u
//             double ue = -blas::DotProd(f,mom) / g_norm;
//             double delta = eps * g_norm / (d-1);
//             double zeta = exp(-delta);
//             for(int i = 0; i < d; ++i)u[i] = (-g[i]/g_norm) *(1-zeta)*(1+zeta + ue * (1-zeta)) + 2*zeta* u[i];

//             //normalize u
//             double u_norm = norm(u, d);
//             for(int i = 0; i < d; ++i)u[i] /= u_norm;
      
//         for(int x=0; x<Nx; x++)
//                 for(int y=0; y<Ny; y++)
//                   for(int mu=0; mu<2; mu++) {
//                 temp = (mom->read(x,y,mu) - f->read(x,y,mu)/f_norm * (sh + ue * (ch - 1)))/ (ch + ue * sh);
//                 mom->write(x,y,mu, temp);

//                   }
      
      
      for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
          for(int mu=0; mu<2; mu++) {
        temp = (mom->read(x,y,mu) - f->read(x,y,mu)/f_norm * (sh + ue * (ch - 1)))/ (ch + ue * sh);
        mom->write(x,y,mu, temp);
    
          }
  }
 
}


void HMC::langevin_noise(field<double> *mom,field<Complex> *gauge){
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  double temp = 0.0;
  int d = Nx * Ny * 2;
  double dtau = gauge->p.tau/gauge->p.n_step;
  double nu = sqrt((exp(2 * dtau / sqrt(d)) - 1.0) / d);
  double langevined_norm2 = 0.;
  double langevined_norm = 0.;
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
        temp = mom->read(x,y,mu) + nu * drand48();
        langevined_norm2 += abs(temp*temp);
        mom->write(x,y,mu, temp);
      }
  langevined_norm = sqrt(langevined_norm2);
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
        temp = mom->read(x,y,mu)/langevined_norm;
        mom->write(x,y,mu, temp);
      }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (f)
void HMC::update_mom_mchmc(field<double> *f, field<double> *mom, double dtau){
  
  int Nx = f->p.Nx;
  int Ny = f->p.Ny;
  int d = Nx * Ny * 2;
  double temp = 0.0;

  int f_norm = sqrt(blas::norm2(f));
  // field<double> *f_normed;
  // blas::copy(f_normed,f);
	// for(int x=0; x<Nx; x++)
	// for(int y=0; y<Ny; y++)
	// for(int mu=0; mu<2; mu++) {
	// f_normed->write(x,y,mu, f_normed->read(x,y,mu)/f_norm);
	// }
  double ue = -1.0 * blas::DotProd(f,mom)/f_norm;
  // double ue = -1.0 * f->dot(mom)/f_norm;
  double sh = sinh(dtau * f_norm / (d-1));
  double ch = cosh(dtau * f_norm / (d-1));
  double th = tanh(dtau * f_norm / (d-1));
  double delta_r = log(ch) + log1p(ue*th);
      
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
          
          
	temp = (mom->read(x,y,mu) + f->read(x,y,mu)/f_norm * (sh + ue * (ch - 1)))/ (ch + ue * sh);
	mom->write(x,y,mu, temp);
          
          
      }
}


//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
// MILC: update_u
void HMC::update_gauge(field<Complex> *gauge, field<double> *mom, double dtau){
  
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

// Optimise this to operate only on a single parity of sites.
// Compute Force = \nabla (\phi^dag Ddag D \phi)
//               = - [(DdagD)^-1 phi^dag] * [(d/dw Ddag) D  + Ddag (d/dw D)] * [(DdagD)^-1 phi]
//
// let (DdagD)^-1 phi = phip
//                 Dp = (d/dw D)
//              Ddagp = (d/dw Ddag)
//
// Force = [ < phip| Ddagp D |phip> - < phip| Ddag Dp |phip >
// Force = [ < phip| Ddagp  |Dphip> - < phipD| Dp |phip >

int HMC::forceD(field<double> *fD, field<Complex> *phi, field<Complex> *gauge)
{
  int cg_iter = 0;
  if(gauge->p.flavours > 0) {

    blas::zero(phip);    
    cg_iter += inv->solve(phip, phi, gauge);
    
    //g3Dphi = g3D * phip
    g3Dpsi(g3Dphi, phip, gauge);

    double fwd_bc = 1.0;
    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
    //#pragma omp parallel for
    for(int x=0; x<Nx; x++) {
      int xp1 = (x+1)%Nx;
      for(int y=0; y<Ny; y++) {
	int yp1 = (y+1)%Ny;
	if(yp1 == 0) fwd_bc = -1.0;
	double temp = 0.0;
	
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
			(fwd_bc*conj(phip->read(x,yp1,0)) * (r*g3Dphi->read(x,y,0) - I*g3Dphi->read(x,y,1)) -
			 fwd_bc*conj(phip->read(x,yp1,1)) * (I*g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
		       -			       
		       (gauge->read(x,y,1) *
			(conj(phip->read(x,y,0)) * fwd_bc * (r*g3Dphi->read(x,yp1,0) + I*g3Dphi->read(x,yp1,1)) +
			 conj(phip->read(x,y,1)) * fwd_bc * (I*g3Dphi->read(x,yp1,0) - r*g3Dphi->read(x,yp1,1))))
		       )
		    );
	
	fD->write(x,y,1,temp);
	fwd_bc = 1.0;
      }
    }
  }
  return cg_iter;
}

// Optimise this to operate only on a single parity of sites.
int HMC::forceMultiD(field<double> *fD, field<Complex> *phi, field<Complex> *gauge)
{
  int n_ferm = force_pfe.inv_pole.size();
  int cg_iter = 0;
  if(gauge->p.flavours == 3) {
    
    std::vector<field<Complex>*> phi_array;
    phi_array.reserve(n_ferm);
    for(int i=0; i<n_ferm; i++) phi_array.push_back(new field<Complex>(gauge->p));    
    blas::zero(fD);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1) phi
    // phip = (D^-1 * g3 * D^-1 g3) phi.
    cg_iter += inv->solveMulti(phi_array, phi, gauge, force_pfe.inv_pole);
    
    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
    //fermion force (D operator)
    field<double> *f = new field<double>(gauge->p);
  
    for(int i=0; i<n_ferm; i++) {
      
      //g3Dphi = g3D * phip
      g3Dpsi(g3Dphi, phi_array[i], gauge);
      blas::zero(f);

      double alpha = force_pfe.inv_res[i];
      double fwd_bc = 1.0;      
      //#pragma omp parallel for
      for(int x=0; x<Nx; x++) {
	int xp1 = (x+1)%Nx;
	for(int y=0; y<Ny; y++) {
	  int yp1 = (y+1)%Ny;
	  if(yp1 == 0) fwd_bc = -1.0;
	  
	  double temp = 0.0;
	  
	  //mu = 0
	  //upper
	  // | r  1 | 
	  // | 1  r |
	  //lower
	  // | r -1 |
	  // | 1 -r |	
	  temp = real(I*((conj(gauge->read(x,y,0)) *
			  (conj(phi_array[i]->read(xp1,y,0)) * (r*g3Dphi->read(x,y,0) +   g3Dphi->read(x,y,1)) -
			   conj(phi_array[i]->read(xp1,y,1)) * (  g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			 -
			 (gauge->read(x,y,0) *
			  (conj(phi_array[i]->read(x,y,0)) * (r*g3Dphi->read(xp1,y,0) -   g3Dphi->read(xp1,y,1)) +
			   conj(phi_array[i]->read(x,y,1)) * (  g3Dphi->read(xp1,y,0) - r*g3Dphi->read(xp1,y,1))))
			 )
		      );
	  f->write(x,y,0,temp*alpha);
	  
	  //mu = 1
	  //upper
	  // | r -i | 
	  // | i  r |
	  //lower
	  // | r  i |
	  // | i -r |
	  temp = real(I*((conj(gauge->read(x,y,1)) *
			  (conj(fwd_bc*phi_array[i]->read(x,yp1,0)) * (r*g3Dphi->read(x,y,0) - I*g3Dphi->read(x,y,1)) -
			   conj(fwd_bc*phi_array[i]->read(x,yp1,1)) * (I*g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			 -			       
			 (gauge->read(x,y,1) *
			  (conj(phi_array[i]->read(x,y,0)) * fwd_bc * (r*g3Dphi->read(x,yp1,0) + I*g3Dphi->read(x,yp1,1)) +
			   conj(phi_array[i]->read(x,y,1)) * fwd_bc * (I*g3Dphi->read(x,yp1,0) - r*g3Dphi->read(x,yp1,1))))
			 )
		      );
	  f->write(x,y,1,temp*alpha);
	  fwd_bc = 1.0;
	}
      }
      // Accumulate force into fD
      // cout << "force size["<<i<<"] = " << blas::norm(f->data) << endl;
      blas::axpy(1.0, f, fD);
    } 
    for(int i=0; i<n_ferm; i++) delete phi_array[i];
    phi_array.resize(0);
    delete f;
  }
  
  return cg_iter;
}

double HMC::measGaugeAction(field<Complex> *gauge) {

  double beta = gauge->p.beta;
  double action = 0.0;
  Complex plaq = 0.0;

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  //#pragma	omp parallel for reduction (+:action)
  for(int x=0; x<Nx; x++) {
    int xp1 = (x+1)%Nx;
    for(int y=0; y<Ny; y++) {
      int yp1 = (y+1)%Ny;
      Complex plaq = (gauge->read(x,y,0) * gauge->read(xp1,y,1) *
		      conj(gauge->read(x,yp1,0)) * conj(gauge->read(x,y,1)));
      
      action += beta*real(1.0 - plaq);      
    }
  }
  return action;
}

double HMC::measMomAction(field<double> *mom) {
  
  double action = 0.0;
  double temp = 0.0;
  int Nx = mom->p.Nx;
  int Ny = mom->p.Ny;
#pragma	omp parallel for reduction (+:action)
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++){
      for(int mu=0; mu<2; mu++){
	double temp = mom->read(x,y,mu);
	action += 0.5 * temp * temp;
      }
    }
  
  return action;
}

//Wilson fermion
double HMC::measFermAction(field<Complex> *gauge, field<Complex> *phi,
			   PFE &pfe, bool rational) {

  int n_ferm = rational ? pfe.inv_pole.size() : 1;
  std::vector<field<Complex>*> phi_tmp;
  phi_tmp.reserve(n_ferm);
  for(int i=0; i<n_ferm; i++) {
    phi_tmp.push_back(new field<Complex>(gauge->p));
    blas::zero(phi_tmp[i]);
  }
  
  Complex action = 0.0;

  if(rational) {
    inv->solveMulti(phi_tmp, phi, gauge, pfe.inv_pole);
    action = pfe.inv_norm * blas::cDotProd(phi, phi);
    for(int i=0; i<n_ferm; i++) action += pfe.inv_res[i] * blas::cDotProd(phi_tmp[i], phi);
  } else {
    inv->solve(phi_tmp[0], phi, gauge);
    action = blas::cDotProd(phi_tmp[0], phi);
  }
  
  for(int i=0; i<n_ferm; i++) delete phi_tmp[i];
  return action.real();
}

//Wilson Action
double HMC::measAction(field<double> *mom, field<Complex> *gauge, std::vector<field<Complex>*> &phi, PFE &pfe) {

  int n_ferm = phi.size();
  double action = 0.0;
  double beta = gauge->p.beta;
  int Ny = gauge->p.Ny;
  int Nx = gauge->p.Nx;
  
  action += measMomAction(mom);
  action += (Nx * Ny)*beta*real(1.0 - measPlaq(gauge));

  if(gauge->p.flavours > 0) action += measFermAction(gauge, phi[0], pfe, false);
  if(gauge->p.flavours == 3) action += measFermAction(gauge, phi[1], pfe, true);
  
  return action;
}


Complex HMC::measPlaq(field<Complex> *gauge) {
  Complex plaq = 0.0;
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  double norm = 1.0/(Nx * Ny);

  //#pragma	omp parallel for reduction(+:plaq)
  for(int x=0; x<Nx; x++) {
    int xp1 = (x+1)%Nx;
    for(int y=0; y<Ny; y++) {
      int yp1 = (y+1)%Ny;
      // Anti-clockwise plaquette, starting at (x,y)
      plaq += (gauge->read(x,y,0) * gauge->read(xp1,y,1) * conj(gauge->read(x,yp1,0)) * conj(gauge->read(x,y,1)));
    }
  }
  return plaq * norm;
}


//----------------------------------------------------------------------------------

HMC::~HMC() {

#ifdef ENABLE_ALG_REMEZ
  delete remez;
#endif
  for(int i=0; i<10; i++) delete guess_stack[i];
  delete phip;
  delete g3Dphi;
  
};
