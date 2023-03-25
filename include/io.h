#pragma once

#include "schwinger2d_internal.h"
#include "blas.h"
#ifdef ENABLE_HDF5
#include "hdf5.h"
#endif

using namespace std;

// Gauge
void writeGauge(field<Complex> *gauge, string name);
void readGauge(field<Complex> *gauge, string name);

// Momentum
void writeMom(field<double> *mom, string name);
void readMom(field<double> *mom, string name);

// PFE
void writePFE(PFE &pfe, string name);
bool readPFE(PFE &pfe, string name);

void hdf5Example();
