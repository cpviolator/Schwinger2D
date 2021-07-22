#pragma once

#include "schwinger2d_internal.h"
#include "blas.h"

#if defined (ENABLE_HDF5)
#include "hdf5.h"
#endif

using namespace std;

void writeGauge(field<Complex> *gauge, string name);
void readGauge(field<Complex> *gauge, string name);
void hdf5Example();
