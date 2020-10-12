#include "hmc.h"
#include "iram.h"

//2D HMC Routines
//---------------------------------------------------------------------
int hmc(field<Complex> *gauge, int iter, double &expdHAve, double &dHAve, int &hmccount) {
  
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
    hmccount++;
    expdHAve += exp(-(H-H_old));
    dHAve += (H-H_old);
  }

  //cout << "H = " << H << "  H old = " << H_old << " delta = " << H - H_old << endl;
  
  // Metropolis accept/reject step
  if (iter >= gauge->p.therm) {    
    if ( drand48() > exp(-(H-H_old)) ) gauge->copy(gauge_old);
    else accept = 1;
  }
  
  return accept;
}

void trajectory(field<double> *mom, field<Complex> *gauge, field<Complex> *phi, int iter){
  
  double dtau = gauge->p.tau/gauge->p.n_step;
  double H = 0.0;
  field<Complex> *guess = new field<Complex>(gauge->p);

  if(gauge->p.deflate) {
    //deflate using phi as source
    //Deflation eigenvectors
    eig_param_t eig_param;
    eig_param.n_ev = gauge->p.n_ev;
    eig_param.n_kr = gauge->p.n_kr;
    eig_param.n_conv = gauge->p.n_conv;
    eig_param.max_restarts = gauge->p.eig_max_restarts;
    eig_param.tol = gauge->p.eig_tol;;
    eig_param.spectrum = 1;
    eig_param.verbose = false;
    
    std::vector<field<Complex> *> kSpace(eig_param.n_conv);
    for(int i=0; i<eig_param.n_conv; i++) kSpace[i] = new field<Complex>(gauge->p);
    //Deflation eigenvalues
    std::vector<Complex> evals(eig_param.n_conv);
    
    guess->copy(phi);
    iram(gauge, kSpace, evals, eig_param);
    deflate(guess, phi, kSpace, evals, eig_param);

    
    
    // Sanity:
    //cout << "Guess norm = " << blas::norm(guess->data) << endl;
    //cout << "phi norm = " << blas::norm(phi->data) << endl;
  }
  
  //gauge force
  field<double> *fU = new field<double>(gauge->p);
  //fermion fermion
  field<double> *fD = new field<double>(gauge->p);
  
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  forceD(fD, gauge, phi, guess);
  update_mom(fU, fD, mom, 0.5*dtau);  
  
  for(int k=1; k<gauge->p.n_step; k++) {
    
    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, dtau);
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge);
    forceD(fD, gauge, phi, guess);
    update_mom(fU, fD, mom, dtau);
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(fU, gauge);
  forceD(fD, gauge, phi, guess);
  update_mom(fU, fD, mom, 0.5*dtau);
  
  //trajectory complete
}

void forceU(field<double> *fU, field<Complex> *gauge) {
  
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
void update_mom(field<double> *fU, field<double> *fD, field<double> *mom, double dtau){
  
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
void update_gauge(field<Complex> *gauge, field<double> *mom, double dtau){
  
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
//----------------------------------------------------------------------------------

