#pragma once

#include "schwinger2d_internal.h"
#include "dirac_op.h"
#include "utils.h"
#include "blas.h"

class IRAM {
  
private:

  // Store an eigendecomposition for comparison methods  
  std::vector<field<Complex>*> kSpace_pre;
  std::vector<Complex> evals_pre;

  // Store a deflation space for CG methods
  std::vector<field<Complex>*> kSpace_defl;
  std::vector<Complex> evals_defl;

  // Store an MG deflation space for CG methods
  std::vector<field<Complex>*> kSpace_mg;
  std::vector<Complex> evals_mg;

  // Object to hold the block orthonormal low mode space
  std::vector<std::vector<Complex>> block_data_ortho;
  // Object to hold the projection coeffiecients of the high modes on the ow space
  std::vector<std::vector<Complex>> block_coeffs;
  
  
  EigParam eig_param; 
  bool verbosity;
  bool use_compressed_space = true;
  int inspection_counter;

  int Nx;
  int Ny;
  
  int n_ev;  
  int n_kr;
  int n_conv;
  int n_deflate;
  int max_restarts;
  double tol;
  bool poly_acc;
  double amax;
  double amin;
  Spectrum spectrum;
  bool iram_verbose;
  Operator op;

  int x_block_size;
  int y_block_size;
  int n_blocks;
  int block_scheme[2];
  int block_size;
  int n_low;
  
public:  
  
  IRAM(EigParam param);

  bool deflationSpaceExists() {return kSpace_defl.size() > 0 ? true : false;};
  
  void OPERATOR(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
  
  void iram(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
	    std::vector<Complex> &evals);

  void computeDeflationSpace(const field<Complex> *gauge);

  void computeMGDeflationSpace(const field<Complex> *gauge);
  
  void deflate(field<Complex> *guess, field<Complex> *phi);
  
  void eigensolveFromUpperHess(MatrixXcd &upperHessEigen, MatrixXcd &Qmat,
			       std::vector<Complex> &evals,
			       std::vector<double> &residua,
			     const double beta, int nKr);

  void qriteration(MatrixXcd &Rmat, MatrixXcd &Qmat, const int nKr, const double tol);
  
  int qrFromUpperHess(MatrixXcd &upperHess, MatrixXcd &Qmat, std::vector<Complex> &evals,
		      std::vector<double> &residua, const double beta, const int nKr,
		      const double tol);
  
  void arnoldiStep(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
		   Eigen::MatrixXcd &upperHessEigen,
		   field<Complex> *r, double &beta, int j);
  
  void reorder(std::vector<field<Complex> *> kSpace, std::vector<Complex> evals, std::vector<double> residua, int nKr, Spectrum spectrum);
  
  void computeEvals(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace, std::vector<double> &residua, std::vector<Complex> &evals, int nEv);
  
  void rotateVecsComplex(std::vector<field<Complex> *> vecs, Eigen::MatrixXcd mat, int num_locked, int iter_keep, int dim);
  
  void zsortc(Spectrum which, int n, std::vector<Complex> &x, std::vector<Complex> &y);
  // Overloaded version of zsortc to deal with real y array.
  void zsortc(Spectrum which, int n, std::vector<Complex> &x, std::vector<double> &y);
  // Overloaded version of zsortc to deal with real x array.
  void zsortc(Spectrum which, int n, std::vector<double> &x, std::vector<Complex> &y);
  // Overloaded version of zsortc to deal with real x and y array.
  void zsortc(Spectrum which, int n, std::vector<double> &x, std::vector<double> &y);
  
  void inspectEvolvedSpectrum(const field<Complex> *gauge, int iter);
  void prepareKrylovSpace(std::vector<field<Complex>*> &kSpace, std::vector<Complex> &evals, const Param p);

  // Read the block data from the (iEig)th vector in kSpace 
  void readVectorToBlock(std::vector<std::vector<Complex>> &block_data);
  
  void blockCompress();
  
  void blockExpand();
  
  ~IRAM() {
    for (int i=0; i<kSpace_pre.size(); i++) delete kSpace_pre[i];
    for (int i=0; i<kSpace_defl.size(); i++) delete kSpace_defl[i];
    for (int i=0; i<kSpace_mg.size(); i++) delete kSpace_mg[i];
    kSpace_pre.resize(0);
    kSpace_defl.resize(0);
    kSpace_mg.resize(0);
    evals_pre.resize(0);    
    evals_defl.resize(0);
    evals_mg.resize(0);
  };
  
};
