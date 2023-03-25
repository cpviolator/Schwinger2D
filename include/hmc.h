#pragma once

#include "schwinger2d_internal.h"
#include "utils.h"
#include "dirac_op.h"
#include "measurements.h"

class HMC {
  
private:

  /// The inverter
  inverterCG *inv;

  // Objects for fermion force and guess tracking
  field<Complex> *phip;
  field<Complex> *g3Dphi;
  std::vector<field<Complex>*> guess_stack;  
  int guess_counter;

#ifdef ENABLE_ALG_REMEZ
  AlgRemez *remez;
#endif
  PFE heatbath_pfe;
  PFE force_pfe;

  // Backup for MCHMC
  field<double> *mom_old;
  
public:

  int hmc_count = 0;
  double exp_dH_ave = 0.0;
  double dH_ave = 0.0;

  double exp_dH = 0.0;
  double dH = 0.0;
  
  HMC(Param param);
  int hmc(field<Complex> *gauge, int iter);
  bool hmc_reversibility(field<Complex> *gauge, int iter);
  void trajectory(field<double> *mom, field<Complex> *gauge, std::vector<field<Complex>*> &phi, int iter);

  int forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi);
  
  void update_mom(field<double> *fU, field<double> *fD, field<double> *mom, double dtau);
  void update_mom(field<double> *f, field<double> *mom, double dtau);
  void update_gauge(field<Complex> *gauge, field<double> *mom, double dtau);
  void forceGradient(field<double> *mom, std::vector<field<Complex>*> phi, field<Complex> *gauge, double one_minus_2lambda_dt, double xi_dtdt);
  void forceGradient(field<double> *mom, field<Complex> *gauge, double one_minus_2lambda_dt, double xi_dtdt);
  void innerFGI(field<double> *mom, field<Complex> *gauge, double tau, int steps);
  void computeFermionForce(field<double>* fD, field<Complex> *gauge, std::vector<field<Complex>*> &phi);
  void computeGaugeForce(field<double>* fU, field<Complex> *gauge);
  
  void langevin_noise(field<double> *mom,field<Complex> *gauge);
    
  double measGaugeAction(field<Complex> *gauge);
  double measMomAction(field<double> *mom);
  double measFermAction(field<Complex> *gauge, field<Complex> *phi, PFE &pfe, bool rational);
  double measAction(field<double> *mom, field<Complex> *gauge, std::vector<field<Complex>*> &phi, PFE &pfe);
  Complex measPlaq(field<Complex> *gauge);
  
  // Optimise this to operate only on a single parity of sites.
  int forceMultiD(field<double> *fD, field<Complex> *phi, field<Complex> *gauge);
  
  ~HMC();
  
};
