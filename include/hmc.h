#pragma once

#include "schwinger2d_internal.h"
#include "utils.h"
#include "dirac_op.h"
#include "iram.h"
#include "measurements.h"

class leapfrogHMC {
  
private:

  // Deflation tracking
  eig_param_t eig_param;
  std::vector<field<Complex>*> kSpace0;
  std::vector<Complex> evals0;

  std::vector<field<Complex>*> kSpace1;
  std::vector<Complex> evals1;

  std::vector<field<Complex>*> kSpace_delta;
  std::vector<Complex> evals_delta;

  std::vector<field<Complex>*> kSpace_prediction;
  std::vector<Complex> evals_prediction;

  // The inverter
  inverterCG *inv;

  // Objects for fermion force and guess tracking
  field<Complex> *phip;
  field<Complex> *g3Dphi;
  std::vector<field<Complex>*> guess_stack;  
  int guess_counter;
  
  
public:

  int hmc_count = 0;
  double exp_dH_ave = 0.0;
  double dH_ave = 0.0;
  
  leapfrogHMC(param_t param);
  int hmc(field<Complex> *gauge, int iter);
  void trajectory(field<double> *mom, field<Complex> *gauge, field<Complex> *phi, int iter);
  void forceU(field<double> *fU, field<Complex> *gauge);
  int forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi,
	     std::vector<field<Complex>*> &kSpace, std::vector<Complex> &evals, int iter);
  void update_mom(field<double> *fU, field<double> *fD, field<double> *mom, double dtau);
  void update_gauge(field<Complex> *gauge, field<double> *mom, double dtau);
  
  ~leapfrogHMC();
  
};
