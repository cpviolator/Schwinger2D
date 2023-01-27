#pragma once

#include "schwinger2d_internal.h"
#include "utils.h"
#include "blas.h"

// Wilson stencil
// D_{W}(n,m) = (m_{0} + 2r)\delta(n,m)
//               - 1/2 Sum [(r-\sigma_{mu}) U_{n,\mu} \delta_{n,m-\hat{\mu}} +
//                          (r+\sigma_{mu}) U^{\dagger}_{m,\mu} \delta_{n,m+\hat{\mu}}]
//
// sigma_1 = | 0  1 |  sigma_2 = | 0 -i | sigma_3 = i*sigma_1*sigma_2 = | 1  0 |
//           | 1  0 |            | i  0 |                               | 0 -1 |
void Dpsi(field<Complex> *psi2, const field<Complex> *psi1, const field<Complex> *gauge);

void g3Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
void Ddagpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
void DdagDpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
void DDdagpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);

void g3psi(field<Complex> *out, const field<Complex> *in);
void g2psi(field<Complex> *out, const field<Complex> *in);
void g1psi(field<Complex> *out, const field<Complex> *in);
