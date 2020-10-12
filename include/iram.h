#pragma once

#include "schwinger2d_internal.h"
#include "dirac_op.h"
#include "utils.h"
#include "blas.h"

void iram(field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
	  std::vector<Complex> &evals, eig_param_t param);

void deflate(field<Complex> *guess, field<Complex> *phi,
	     std::vector<field<Complex> *> kSpace, std::vector<Complex> &evals,
	     eig_param_t param);

void eigensolveFromUpperHess(MatrixXcd &upperHessEigen, MatrixXcd &Qmat,
			    std::vector<Complex> &evals,
			    std::vector<double> &residua,
			     const double beta, int nKr);

void qriteration(MatrixXcd &Rmat, MatrixXcd &Qmat, const int nKr, const double tol);

int qrFromUpperHess(MatrixXcd &upperHess, MatrixXcd &Qmat, std::vector<Complex> &evals,
		    std::vector<double> &residua, const double beta, const int nKr,
		    const double tol);

void arnoldiStep(field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
		 Eigen::MatrixXcd &upperHessEigen,
		 field<Complex> *r, double &beta, int j);

void reorder(std::vector<field<Complex> *> kSpace, std::vector<Complex> evals, std::vector<double> residua, int nKr, int spectrum);

void computeEvals(field<Complex> *gauge, std::vector<field<Complex> *> kSpace, std::vector<double> &residua, std::vector<Complex> &evals, int nEv);

void rotateVecsComplex(std::vector<field<Complex> *> vecs, Eigen::MatrixXcd mat, int num_locked, int iter_keep, int dim);

void zsortc(int which, int n, std::vector<Complex> &x, std::vector<Complex> &y);
// Overloaded version of zsortc to deal with real y array.
void zsortc(int which, int n, std::vector<Complex> &x, std::vector<double> &y);
// Overloaded version of zsortc to deal with real x array.
void zsortc(int which, int n, std::vector<double> &x, std::vector<Complex> &y);
// Overloaded version of zsortc to deal with real x and y array.
void zsortc(int which, int n, std::vector<double> &x, std::vector<double> &y);
