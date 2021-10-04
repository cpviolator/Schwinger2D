#include "dirac_op.h"

// X = (x,y)
// a,b are spin indices
// mu is gamma/direction index
// \sigma_1 = 0  1
//            1  0 
//
// \sigma_2 = 0 -i
//          = i  0
// 
// \sigma_3 = 1  0
//            0 -1

// D_{wilson} = (m0 + 2) * I_{X}
//              - 1/2 \Sum_{mu}   (I_{a,b} - \sigma_{\mu})U(X)_{\mu} \delta_{X+\mu}
//                              + (I_{a,b} + \sigma_{\mu})U^*(X-mu)_{\mu}) \delta_{X-mu}
//
void Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge) {
  
  double m0 = gauge->p.m;
  double  r = 1.0;
  double constant = (m0 + 2*r);
  double fwd_bc = 1.0;
  double bwd_bc = 1.0;  
  
  //Sum over 0,1 directions at each lattice site.
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;

  //#pragma omp parallel for collapse(2)
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {

      int xp1 = (x+1)%Nx;
      int xm1 = (x-1+Nx)%Nx;    
      int yp1 = (y+1)%Ny;
      int ym1 = (y-1+Ny)%Ny;

      if(yp1 == 0) fwd_bc = -1.0;
      if(y   == 0) bwd_bc = -1.0;  
      
      
      Complex tmp = 0.0;
      
      //upper
      //(m0 + 2) * I_{X}
      tmp = constant * in->read(x,y,0)

	// mu = 1
	- 0.5*(gauge->read(x,y,0)         * (in->read(xp1,y,0) - in->read(xp1,y,1)) +
	       conj(gauge->read(xm1,y,0)) * (in->read(xm1,y,0) + in->read(xm1,y,1)))
	
	// mu = 2
	- 0.5*(gauge->read(x,y,1)         * fwd_bc*(in->read(x,yp1,0) + I*in->read(x,yp1,1)) +
	       conj(gauge->read(x,ym1,1)) * bwd_bc*(in->read(x,ym1,0) - I*in->read(x,ym1,1)));
      
      out->write(x,y,0, tmp);
      
      //lower
      //(m0 + 2) * I_{X}
      tmp = constant * in->read(x,y,1) 

	// mu = 1
	- 0.5*(gauge->read(x,y,0)         * (in->read(xp1,y,1) - in->read(xp1,y,0)) +
	       conj(gauge->read(xm1,y,0)) * (in->read(xm1,y,1) + in->read(xm1,y,0)))
	
	// mu = 2
	- 0.5*(gauge->read(x,y,1)         * fwd_bc*(in->read(x,yp1,1) - I*in->read(x,yp1,0)) +
	       conj(gauge->read(x,ym1,1)) * bwd_bc*(in->read(x,ym1,1) + I*in->read(x,ym1,0)));
      
      out->write(x,y,1, tmp);

      fwd_bc = 1.0;
      bwd_bc = 1.0;  
    }
  }
}

void g3psi(field<Complex> *out, const field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
  //#pragma omp parallel for 
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0,  in->read(x,y,0));
      out->write(x,y,1, -in->read(x,y,1));
    }
}

void g2psi(field<Complex> *out, const field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
  //#pragma omp parallel for 
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0, -I*in->read(x,y,1));
      out->write(x,y,1,  I*in->read(x,y,0));
    }
}

void g1psi(field<Complex> *out, field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
  //#pragma omp parallel for 
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0, in->read(x,y,1));
      out->write(x,y,1, in->read(x,y,0));
    }
}

void g3Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){
  
  field<Complex> *temp = new field<Complex>(in->p);  
  Dpsi(temp, in, gauge);
  g3psi(out, temp);
  delete temp;
}

void Ddagpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){
  
  field<Complex> *temp = new field<Complex>(in->p);  
  g3psi(out, in);
  Dpsi(temp, out, gauge);
  g3psi(out, temp);
  delete temp;
}


void DdagDpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge) {
  
  field<Complex> *temp = new field<Complex>(in->p);
  Dpsi(temp, in, gauge);
  g3psi(out, temp);
  Dpsi(temp, out, gauge);
  g3psi(out, temp);
  delete temp;
}
