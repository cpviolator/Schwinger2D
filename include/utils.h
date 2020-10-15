#pragma once

#include "schwinger2d_internal.h"
#include "blas.h"
#include "hdf5.h"

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
  bool verbosity = true;
  
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
  int n_conv = 16;
  int n_deflate = 16;
  double eig_tol = 1e-6;
  int eig_max_restarts = 10000;
  bool poly_acc = false;
  double amax = -1.0;
  double amin = -1.0;
  int poly_deg = 0;

  //Eigensolver compression
  int block_scheme[2];
  int n_low = 4;
  
  //Measurements
  bool meas_pl = false; //Polyakov loops
  bool meas_wl = false; //Wilson loop and Creutz ratios
  bool meas_pc = false; //Pion
  bool meas_vt = false; //Vacuum trace

  //Wilson loop and Polyakov loop max size.
  int loop_max = 16;
  
} param_t;

typedef struct {

  int n_ev = 0;
  int n_kr = 0;
  int n_conv = 0;
  int n_deflate = 0;
  int max_restarts = 0;
  double tol = 0.0;
  int spectrum = 0;
  bool verbose = false;

  int block_size = 0;
  
} eig_param_t;


template<typename T> class field {
  
public:
  
  ~field(){
    data.resize(0);
  }
  
  std::vector<T> data;
  param_t p;
  
  field(std::vector<T> &data, param_t p) : data(data), p(p) { }    
  field(param_t p) : p(p)
  {
    data.resize(p.Nx * p.Ny * 2);
    for(unsigned int i=0; i<data.size(); i++) data[i] = 0.0;
  }
  
  T read(int x, int y, int mu) {
    return data[2*(x + p.Nx * y) + mu];
  }
  
  void write(int x, int y, int mu, const T elem) {
    data[2*(x + p.Nx * y) + mu] = elem;
  }
  
  void copy(field<T> *in){
    blas::copy(data, in->data);
  }

  unsigned int size() { return data.size(); }
  
  void print() {
    for(int x=0; x<p.Nx; x++) {
      for(int y=0; y<p.Ny; y++) {
	for(int mu=0; mu<2; mu++) {      
	  cout << "elem("<<x<<","<<y<<":" << mu << ") = " << data[2*(x + p.Nx * y) + mu] << endl;
	}
      }
    }    
  }
  
};

void printParams(param_t p);
void constructName(string &name, param_t p);
void writeGauge(field<Complex> *gauge, string name);
void readGauge(field<Complex> *gauge, string name);

void gaussStart(field<Complex> *gauge);
void coldStart(field<Complex> *gauge);

// Normalized gaussian exp(-phi*phi/2) and  <phi|phi> = 1
void gaussReal(field<double> *field);

//normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
void gaussComplex(field<Complex> *field);

//APE smearing: project back on U(1)
// staple x is 0th, y is 1st.
void smearLink(field<Complex> *smeared, field<Complex> *gauge);

void measBlockColinearity(std::vector<field<Complex> *> kSpace, int blockScheme[2], int nLow);

void blockCompress(std::vector<field<Complex> *> &kSpace,
		   std::vector<std::vector<Complex>> &block_data_ortho,
		   std::vector<std::vector<Complex>> &block_coef,
		   int blockScheme[2], int n_low, int n_conv);

void blockExpand(std::vector<field<Complex> *> &kSpace,
		 std::vector<std::vector<Complex>> &block_data_ortho,
		 std::vector<std::vector<Complex>> &block_coef,
		 int blockScheme[2], int n_low, int n_conv);

void hdf5Example();
