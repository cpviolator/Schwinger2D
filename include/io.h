#pragma once

#include "schwinger2d_internal.h"
#include "blas.h"
#include "hdf5.h"

using namespace std;

void writeGauge(field<Complex> *gauge, string name);
void readGauge(field<Complex> *gauge, string name);
void hdf5Example();
