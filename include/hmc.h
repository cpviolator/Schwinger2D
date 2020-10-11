#pragma once

#include "schwinger2d_internal.h"
#include "utils.h"
#include "dirac_op.h"
#include "measurements.h"

void trajectory(field<double> *mom, field<Complex> *gauge, field<Complex> *phi, int iter);
int hmc(field<Complex> *gauge, int iter, double &expdHave, double &dHave, int &hmccount);
void forceU(field<double> *fU, field<Complex> *gauge);
void update_mom(field<double> *fU, field<double> *fD, field<double> *mom, double dtau);
void update_gauge(field<Complex> *gauge, field<double> *mom, double dtau);
