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

  // Compute a the low spectrum
  if(iter >= 2*gauge->p.therm && inspectrum_bool) inspectrum(gauge, iter);

  // Construct objects for an eigensolver
  //-------------------------------------
  eig_param_t eig_param;
  std::vector<field<Complex>*> kSpace;
  std::vector<Complex> evals;	
  prepareKrylovSpace(kSpace, evals, eig_param, gauge->p);
  
  // Compute a deflation space using IRAM
  //if(iter >= 2*gauge->p.therm && gauge->p.deflate) iram(gauge, kSpace, evals, eig_param);  
  
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

    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge);
    ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);    
    update_mom(fU, fD, mom, dtau);
  }

  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  ave_iter += forceD(fD, phi, gauge, kSpace, evals, iter);
  update_mom(fU, fD, mom, 0.5*dtau);
  
  // HMC trajectory complete
  //----------------------------------------------------------
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

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
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
    //for(int i=0; i<10; i++) phip->print(i);
    
    //g3Dphi = g3D * phip
    g3Dpsi(g3Dphi, phip, gauge);
    
    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
#pragma omp parallel for
    for(int x=0; x<Nx; x++) {
      int xp1 = (x+1)%Nx;
      for(int y=0; y<Ny; y++) {
	int yp1 = (y+1)%Ny;

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
