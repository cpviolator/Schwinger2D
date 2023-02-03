#pragma once

//C++ standard
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <sys/time.h>

//Externals
#ifdef ENABLE_ALG_REMEZ
#include "alg_remez.h"
#endif

#include <Eigen/Eigenvalues>
using namespace std;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;

//Convenient typedef for complex double
typedef complex<double> Complex;

#define PI 3.141592653589793
#define TWO_PI 6.283185307179586
#define I Complex(0,1.0)
#define cUnit Complex(1.0,0)

enum Integrator { LEAPFROG = 0, FGI = 1};
enum Operator { M = 0, Mdag = 1, MMdag = 2, MdagM = 3};
// (L)argest (S)mallest (R)eal (I)maginary
enum Spectrum { LM = 0, SM = 1, LR = 2, SR = 3, LI = 4, SI = 6};

using namespace std;

typedef struct EigParam {

  int Nx = 32;
  int Ny = 32;
  
  int n_ev = 16;
  int n_kr = 32;
  int n_conv = 16;
  int n_deflate = 0;
  int max_restarts = 100;
  double tol = 1e-9;
  bool poly_acc = false;
  int poly_deg = 32;
  double amax = -1.0;
  double amin = -1.0;
  Spectrum spectrum = SR;
  bool verbosity = false;
  bool iram_verbose = false;
  Operator op = MdagM;  
  int block_scheme[2];
  int n_low = 16;
  
} EigParam_t;

typedef struct PFE {
  
  double norm = 0;
  std::vector<double> res;
  std::vector<double> pole;

  double inv_norm = 0;
  std::vector<double> inv_res;
  std::vector<double> inv_pole;
  
} PFE_t;


class Param{
  
public:

  int seed = 1234;
  bool verbosity = true;
  int current_hmc_iter = 0;
  
  // Physics
  int Nx = 32;
  int Ny = 32;
  double beta = 3.0;
  double m = 0.1;
  double m_heavy = 0.5;
  int flavours = 2;
  
  // HMC
  int n_step = 4;
  int inner_step = 1;
  double tau = 1.0;
  int iter_hmc = 1000;
  int therm = 250;
  int skip = 5;
  int chkpt = 100;
  int checkpoint_start = 0;
  int pfe_degree = 15;
  int pfe_prec = 50;  
  int max_iter_cg = 1000;
  double tol_cg = 1e-9;
  bool cg_verbosity = false;
  Integrator integrator = FGI;
  int reverse = 100;
      
  // Eigensolver params
  EigParam eig_param;
  bool inspect_spectrum = true;
  bool deflate = false;

  //Smearing
  double alpha = 0.5;
  int smear_iter = 1;
  
  //Measurements
  bool meas_pl = false; //Polyakov loops
  bool meas_wl = false; //Wilson loop and Creutz ratios
  bool meas_pc = false; //Pion
  bool meas_vt = false; //Vacuum trace

  //Wilson loop and Polyakov loop max size.
  int loop_max = std::min(Ny/2, Nx/2);
  
  void usage(char **argv);
  void print();
  int init(int argc, char **argv, int *idx);
  
};

