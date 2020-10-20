#include "dirac_op.h"

void Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge) {
  
  double m0 = gauge->p.m;
  double  r = 1.0;
  double constant = (2*r + m0);

  // Sanity checks
  //cout << "Norm in = " << blas::norm(in->data);
  
  //Sum over 0,1 directions.
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
#pragma omp parallel for 
  for(int x=0; x<Nx; x++) {
    int xp1 = (x+1)%Nx;
    int xm1 = (x-1+Nx)%Nx;    
    for(int y=0; y<Ny; y++) {
      int yp1 = (y+1)%Ny;
      int ym1 = (y-1+Ny)%Ny;            
      //upper
      Complex tmp = 0.0;
      tmp = constant * in->read(x,y,0) -
	
	0.5*(     gauge->read(x,y,0)    * (r*in->read(xp1,y,0) - in->read(xp1,y,1)) +
	     conj(gauge->read(xm1,y,0)) * (r*in->read(xm1,y,0) + in->read(xm1,y,1)) +
		  
		  gauge->read(x,y,1)    * (r*in->read(x,yp1,0) + I*in->read(x,yp1,1)) +
	     conj(gauge->read(x,ym1,1)) * (r*in->read(x,ym1,0) - I*in->read(x,ym1,1)));
      out->write(x,y,0, tmp);
      
      //lower
      tmp = constant * in->read(x,y,1) -
	
	0.5*(     gauge->read(x,y,0)    * (-in->read(xp1,y,0) + r*in->read(xp1,y,1)) -
	     conj(gauge->read(xm1,y,0)) * (-in->read(xm1,y,0) - r*in->read(xm1,y,1)) +

		  gauge->read(x,y,1)    * (-I*in->read(x,yp1,0) + r*in->read(x,yp1,1)) -
	     conj(gauge->read(x,ym1,1)) * (-I*in->read(x,ym1,0) - r*in->read(x,ym1,1)));
      out->write(x,y,1, tmp);
    }
  }
  // Sanity
  //cout << "Norm out = " << blas::norm(out->data) << endl;
  //exit(0);
}

void g3psi(field<Complex> *out, const field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
#pragma omp parallel for 
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0,  in->read(x,y,0));
      out->write(x,y,1, -in->read(x,y,1));
    }
}

void g2psi(field<Complex> *out, const field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
#pragma omp parallel for 
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0, -I*in->read(x,y,1));
      out->write(x,y,1,  I*in->read(x,y,0));
    }
}

void g1psi(field<Complex> *out, field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
#pragma omp parallel for 
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
