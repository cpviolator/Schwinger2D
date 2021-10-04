#pragma once

#include "schwinger2d_internal.h"
#include "blas.h"
#include "hdf5.h"

using namespace std;

// Gauge
void writeGauge(field<Complex> *gauge, string name);
void readGauge(field<Complex> *gauge, string name);

// PFE
void writePFE(PFE &pfe, string name);
bool readPFE(PFE &pfe, string name);

void hdf5Example();
