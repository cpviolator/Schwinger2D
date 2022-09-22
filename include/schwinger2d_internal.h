#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <sys/time.h>

#include <Eigen/Eigenvalues>
using namespace std;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;

typedef complex<double> Complex;

#define PI 3.141592653589793
#define TWO_PI 6.283185307179586
#define I Complex(0,1.0)
#define cUnit Complex(1.0,0)

using namespace std;

typedef struct param {
  
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
  int inspect_spectrum = false;

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

typedef struct eig_param {

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
