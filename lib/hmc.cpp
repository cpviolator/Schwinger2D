#include "hmc.h"
#include "inverters.h"
#include "iram.h"

//2D HMC Routines
//---------------------------------------------------------------------

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
  
  if (iter >= gauge->p.therm) H_old = measAction(mom, gauge, chi, false);
  trajectory(mom, gauge, phi, iter);
  if (iter >= gauge->p.therm) H = measAction(mom, gauge, phi, true);
  
  if (iter >= 2*gauge->p.therm) {      
    hmc_count++;
    exp_dH_ave += exp(-(H-H_old));
    dH_ave += (H-H_old);
  }

  //cout << "H = " << H << "  H old = " << H_old << " delta = " << H - H_old << endl;
  
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
  
  //gauge force (U field)
  field<double> *fU = new field<double>(gauge->p);
  //fermion force (D operator)
  field<double> *fD = new field<double>(gauge->p);

  // Compute a the low spectrum
  if(iter >= 2*gauge->p.therm && inspectrum_bool) inspectrum(gauge, iter);

  // Start HMC trajectory
  //----------------------------------------------------------
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  forceD(fD, gauge, phi);
  update_mom(fU, fD, mom, 0.5*dtau);  
  
  for(int k=1; k<gauge->p.n_step; k++) {
    
    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, dtau);
    
    // Compute the low spectrum
    if(iter >= 2*gauge->p.therm && inspectrum_bool) inspectrum(gauge, iter);
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge);
    forceD(fD, gauge, phi);
    
    update_mom(fU, fD, mom, dtau);
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau);
  
  // Compute the low spectrum
  if(iter >= 2*gauge->p.therm && inspectrum_bool) inspectrum(gauge, iter);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  forceD(fD, gauge, phi);
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
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {

      xp1 = (x+1)%Nx;
      xm1 = (x-1+Nx)%Nx;
      yp1 = (y+1)%Ny;
      ym1 = (y-1+Ny)%Ny;
      
      plaq0 = gauge->read(x,y,0)*gauge->read(xp1,y,1)*conj(gauge->read(x,yp1,0))*conj(gauge->read(x,y,1));
      plaq =  gauge->read(x,ym1,0)*gauge->read(xp1,ym1,1)*conj(gauge->read(x,y,0))*conj(gauge->read(x,ym1,1));
      temp = beta*(imag(plaq0) - imag(plaq));
      fU->write(x,y,0, temp);
      
      plaq =  gauge->read(x,y,1)*conj(gauge->read(xm1,yp1,0))*conj(gauge->read(xm1,y,1))*gauge->read(xm1,y,0);
      temp = beta*(imag(plaq) - imag(plaq0));
      fU->write(x,y,1, temp);

      //This plaquette was aleady computed. We want the conjugate.
      //fU->read(x,y,1) -= p.beta*imag(plaq0);      
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

/*
field<Complex> *guess_exp;
field<Complex> *sol1;
field<Complex> *sol2;
field<Complex> *delta_sol;

bool init = false;
int counter = 0;
//bool spline = true;
bool spline = false;

// Krylov space
std::vector<field<Complex>*> kSpace(eig_param.n_conv);
for(int i=0; i<eig_param.n_conv; i++) kSpace[i] = new field<Complex>(gauge->p);
// eigenvalues
std::vector<Complex> evals(eig_param.n_conv);
*/

// let dD = (d/dtheta D) :
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
// We should optimise this to operate only on a single parity of sites.
//
// input 
void leapfrogHMC::forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi){
  
  if(gauge->p.dynamic == true) {

    blas::zero(fD->data);
    
    //phip = (D^dagD)^-1 * phi
    field<Complex> *phip = new field<Complex>(gauge->p);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1)
    //  phi = (D^-1 * g3 * D^-1 g3) phi.
    field<Complex> *guess = new field<Complex>(gauge->p);
    gaussComplex(guess);
    blas::ax(10.0, guess->data);
    Ainvpsi(phip, phi, guess, gauge);
        
    //g3Dphi = g3D * phip
    field<Complex> *g3Dphi = new field<Complex>(gauge->p); 
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
    if(init == false && spline) {
      guess_exp = new field<Complex>(gauge->p);
      sol1 = new field<Complex>(gauge->p);
      sol2 = new field<Complex>(gauge->p);
      delta_sol = new field<Complex>(gauge->p);
      init = true;
    }
    
    if(counter == 0 && spline) {
      sol1->copy(x);
    }
    
    if(counter > 0 && spline) {
      sol2->copy(x);
      blas::zero(delta_sol->data);
      // delta_sol = sol2 - sol1 
      blas::axpy( 1.0, sol2->data, delta_sol->data);
      blas::axpy(-1.0, sol1->data, delta_sol->data);
      // sol1 updated
      sol1->copy(sol2);
      // guess updated
      guess_exp->copy(sol2);
      blas::axpy(1.0, delta_sol->data, guess_exp->data);
    }
    
    counter++;
    */
    
    delete phip;
    delete g3Dphi;    
  }
}
//----------------------------------------------------------------------------------

