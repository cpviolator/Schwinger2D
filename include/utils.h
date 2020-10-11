#pragma once

#include "schwinger2d_internal.h"
//#include "measurements.h"
#include "blas.h"

using namespace std;

typedef struct{
  
  //HMC
  int n_step = 25;
  double tau = 1.0;
  int iter_hmc = 1000;
  int therm = 50;
  int skip = 25;
  int chkpt = 100;
  int checkpoint_start = 0;
  int max_iter_cg = 1000;
  double eps = 1e-6;

  int seed = 1234;
  
  //physics
  int Nx = 32;
  int Ny = 32;
  double beta = 3.0;
  double m = -0.06;
  bool dynamic = true;
  bool deflate = true;
  
  //Smearing
  double alpha = 0.5;
  int smear_iter = 1;
  
  //Eigensolver params
  int n_ev = 16;
  int n_kr = 64;
  double eig_tol = 1e-6;
  int eig_max_restarts = 10000;
  bool poly_acc = false;
  double amax = -1.0;
  double amin = -1.0;
  int poly_deg = 0;

  //Measurements
  bool meas_pl = false; //Polyakov loops
  bool meas_wl = false; //Wilson loop and Creutz ratios
  bool meas_pc = false; //Pion
  bool meas_vt = false; //Vacuum trace

  //Wilson loop and Polyakov loop max size.
  int loop_max = 16;
  
} param_t;

template<class T> class field {
  
public:

  ~field();
  
  std::vector<T> data;
  param_t p;
  
  field(std::vector<T> &field, param_t p);
  field(param_t p);  
  T read(int x, int y, int mu);
  void write(int x, int y, int mu, const T elem);
  void copy(field<T> *in);
  void print();
  
};

void printParams(param_t p);
void constructName(string &name, param_t p);
void writeGauge(field<Complex> *gauge, string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  //outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge).real() << endl;

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	outPutFile << setprecision(12) <<  setw(20) << arg(gauge->read(x,y,mu)) << endl;
  
  outPutFile.close();
  return;  
}

void readGauge(field<Complex> *gauge, string name){

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  //Header check
  //getline(inPutFile, val);
  //double plaq_real_header = stod(val);
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	getline(inPutFile, val);
	gauge->write(x, y, mu, polar(1.0, stod(val)));
      }
    }
  }
  
  //double plaq_real_measured = measPlaq(gauge).real();
  //double err = fabs(1.0 - plaq_real_header/plaq_real_measured);
  //if(abs(err) > 1e-12) {
  //cout << "Gauge read fail!" << endl;
  //cout << setprecision(16) << setw(20) << "Plaqette on file  = " << plaq_real_header << endl;
  //cout << setprecision(16) << setw(20) << "Plaqette measured = " << plaq_real_measured << endl;   
  //exit(0);
  //}    
  return;
}

void gaussStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	gauge->write(x, y, mu, polar(1.0,drand48()));  
}  

void coldStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	gauge->write(x, y, mu, Complex(1.0,0.0));    
}  

// Normalized gaussian exp(-phi*phi/2) and  <phi|phi> = 1
void gaussReal(field<double> *field) {
  
  double r, theta, sum;
  int Nx = field->p.Nx;
  int Ny = field->p.Ny;
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	
	r = sqrt(-2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field->write(x,y,mu, r*cos(theta));
      }
    }
  }  
}

//normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
void gaussComplex(field<Complex> *field) {

  double r1, theta1, r2, theta2, sum;
  double inv_sqrt2 = 1.0/sqrt(2);

  int Nx = field->p.Nx;
  int Ny = field->p.Ny;
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	
	r1 = sqrt(-2.0*log(drand48()));
	theta1 = TWO_PI*(drand48());
	r2 = sqrt(-2.0*log(drand48()));
	theta2 = TWO_PI*(drand48());
	
	field->write(x,y,mu, Complex(r1*cos(theta1),r2*sin(theta2))*inv_sqrt2);
      }
    }
  }
}

//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(field<Complex> *smeared, field<Complex> *gauge){

  double alpha = gauge->p.alpha;
  int iter = gauge->p.smear_iter;
  Complex tmp = 0;
  
  field<Complex> *smeared_tmp = new field<Complex>(gauge->p);
  smeared->copy(gauge);
  smeared_tmp->copy(smeared);
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int i=0; i<iter; i++) {    
    for(int x=0; x<Nx; x++) {
      for(int y=0; y<Ny; y++) {
	
	//|->-|   |   |
	//^   v + v   ^
	//|   |   |->-|
	tmp = alpha * (smeared->read(x,y,1) * smeared->read(x,y+1,0) * conj(smeared->read(x+1,y,1)));
	
	tmp += alpha * (conj(smeared->read(x,y-1,1)) * smeared->read(x,y-1,0) * smeared->read(x+1,y-1,1));
			
	smeared_tmp->write(x,y,0, tmp);
				
	//|->-    -<-|
	//^    +     ^
	//|-<-    ->-|
	tmp = alpha * (smeared->read(x,y,0) * smeared->read(x+1,y,1) * conj(smeared->read(x,y+1,0)));
	
	tmp += alpha * (conj(smeared->read(x-1,y,0)) * smeared->read(x-1,y,1) * smeared->read(x-1, y+1,01));
	
	smeared_tmp->write(x,y,1, tmp);	
      }
    }
    
    //Project back to U(1)
    for(int x=0; x<Nx; x++)
      for(int y=0; y<Ny; y++)
	for(int mu=0; mu<2; mu++)
	  smeared->write(x,y,mu, polar(1.0,arg(smeared_tmp->read(x,y,mu))));
  }
}
