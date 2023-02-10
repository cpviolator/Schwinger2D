#include "eig.h"

#ifdef USE_FEAST
#include "feast.h"
#include "feast_sparse.h"
#endif

int int_round(double d)
{
  return (int)floor(d + 0.5);
}

Eig::Eig(EigParam eig_param_in) {
  
  eig_param = eig_param_in;
  verbosity = eig_param.verbosity;
  use_compressed_space = eig_param.use_comp_space;
    
  inspection_counter = 0;
  
  Nx = eig_param.Nx;
  Ny = eig_param.Ny;
  
  n_ev = eig_param.n_ev;
  n_kr = eig_param.n_kr;
  n_conv = eig_param.n_conv;
  n_deflate = eig_param.n_deflate;
  max_restarts = eig_param.max_restarts;
  tol = eig_param.tol;
  spectrum = eig_param.spectrum;
  eig_verbose = eig_param.iram_verbose;
  op = eig_param.op;
  block_scheme[0] = eig_param.block_scheme[0];
  block_scheme[1] = eig_param.block_scheme[1];
  n_low = eig_param.n_low;

  x_block_size = Nx/block_scheme[0];
  y_block_size = Ny/block_scheme[1];
  n_blocks = block_scheme[0]*block_scheme[1];
  block_size = 2 * x_block_size * y_block_size;

  feast_Ncontour = eig_param.feast_Ncontour;
  feast_M0 = eig_param.feast_M0;
  feast_Emax = eig_param.feast_Emax;
  
  if(Nx%block_scheme[0] != 0) {
    cout << "Eig: Error: x block_scheme = " << block_scheme[0] << " does not divide Nx = " << Nx << endl;
    exit(0);
  }
  if(Ny%block_scheme[1] != 0) {
    cout << "Eig: Error: y block_scheme = " << block_scheme[1] << " does not divide Ny = " << Ny << endl;
    exit(0);
  }

  // n_blocks * n_conv * (block_size + l_low) complex elems
  for(int i=0; i<n_blocks; i++) {
    block_data_ortho.push_back(std::vector<Complex> (n_low * block_size, 0.0));
    block_coeffs.push_back(std::vector<Complex> (n_low * n_conv, 0.0));
  }
}

void Eig::OPERATOR(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){
  switch(op) {
  case M: Dpsi(out, in, gauge); break;
  case Mdag: Ddagpsi(out, in, gauge); break;
  case MdagM: DdagDpsi(out, in, gauge); break;
  case MMdag: DDdagpsi(out, in, gauge); break;
  default: cout << "Eig: Undefined operator type requested: " << op << endl;
    exit(0);
  }
}

void Eig::deflate(field<Complex> *guess, field<Complex> *phi) {

  if(kSpace_defl.size() == 0) {
    cout << "Eig Error: No deflation space evecs exist" << endl;
    exit(0);
  }
  if(evals_defl.size() == 0) {
    cout << "Eig: Error: No deflation space evals exist" << endl;
    exit(0);
  }
  
  blas::zero(guess);
  Complex scalar;
  
  // Deflate each converged eigenpair from the guess  
  // guess_defl = (v * lambda^-1 * v^dag) * guess
  if(verbosity) cout << "Eig: Deflating " << n_deflate << " eigenpairs." << endl;
  for(int i=0; i<n_deflate; i++) {
    
    if(use_compressed_space) { 
      // Compute scalar part: s = (lambda)^-1 * (v^dag * phi)
      scalar = blas::cDotProd(kSpace_mg[i], phi);
      scalar /= real(evals_mg[i]);
      // Accumulate in guess defl_guess: defl_guess = defl_guess + v * s
      blas::axpy(scalar, kSpace_mg[i], guess);
    }
    else {
      // Compute scalar part: s = (lambda)^-1 * (v^dag * phi)
      scalar = blas::cDotProd(kSpace_defl[i], phi);
      scalar /= real(evals_defl[i]);
      // Accumulate in guess defl_guess: defl_guess = defl_guess + v * s
      blas::axpy(scalar, kSpace_defl[i], guess);
    }
  }
}

void Eig::zsortc(Spectrum which, int n, std::vector<Complex> &x, std::vector<Complex> &y) {
  
  std::vector<std::pair<Complex, Complex>> array(n);
  for(int i=0; i<n; i++) array[i] = std::make_pair(x[i], y[i]);
  
  switch(which) {
  case LM: std::sort(array.begin(), array.begin()+n,
		    [] (const pair<Complex,Complex> &a,
			const pair<Complex,Complex> &b) {
		      return (abs(a.first) < abs(b.first)); } );
    break;
  case SM: std::sort(array.begin(), array.begin()+n,
		    [] (const pair<Complex,Complex> &a,
			const pair<Complex,Complex> &b) {
		      return (abs(a.first) > abs(b.first)); } );
    break;
  case LR: std::sort(array.begin(), array.begin()+n,
		    [] (const pair<Complex,Complex> &a,
			const pair<Complex,Complex> &b) {
		      return (a.first).real() < (b.first).real(); } );
    break;
  case SR: std::sort(array.begin(), array.begin()+n,
		    [] (const pair<Complex,Complex> &a,
			const pair<Complex,Complex> &b) {
		      return (a.first).real() > (b.first).real(); } );
    break;
  case LI: std::sort(array.begin(), array.begin()+n,
		    [] (const pair<Complex,Complex> &a,
			const pair<Complex,Complex> &b) {
		      return (a.first).imag() < (b.first).imag(); } );
    break;
  case SI: std::sort(array.begin(), array.begin()+n,
		    [] (const pair<Complex,Complex> &a,
			const pair<Complex,Complex> &b) {
		      return (a.first).imag() > (b.first).imag(); } );
    break;
  default: cout << "Eig: Undefined sort" << endl;
  }

  // Repopulate x and y arrays with sorted elements
  for(int i=0; i<n; i++) {
    x[i] = array[i].first;
    y[i] = array[i].second;
  }
}

// Overloaded version of zsortc to deal with real y array.
void Eig::zsortc(Spectrum which, int n, std::vector<Complex> &x, std::vector<double> &y) {

  std::vector<Complex> y_tmp(n,0.0);
  for(int i=0; i<n; i++) y_tmp[i].real(y[i]);
  zsortc(which, n, x, y_tmp);
  for(int i=0; i<n; i++) y[i] = y_tmp[i].real();
}

// Overloaded version of zsortc to deal with real x array.
void Eig::zsortc(Spectrum which, int n, std::vector<double> &x, std::vector<Complex> &y) {

  std::vector<Complex> x_tmp(n,0.0);
  for(int i=0; i<n; i++) x_tmp[i].real(x[i]);
  zsortc(which, n, x_tmp, y);
  for(int i=0; i<n; i++) x[i] = x_tmp[i].real();
}

// Overloaded version of zsortc to deal with real x and y array.
void Eig::zsortc(Spectrum which, int n, std::vector<double> &x, std::vector<double> &y) {

  std::vector<Complex> x_tmp(n,0.0);
  std::vector<Complex> y_tmp(n,0.0);
  for(int i=0; i<n; i++) {
    x_tmp[i].real(x[i]);
    y_tmp[i].real(y[i]);
  }
  zsortc(which, n, x_tmp, y_tmp);
  for(int i=0; i<n; i++) {
    x[i] = x_tmp[i].real();
    y[i] = y_tmp[i].real();
  }
}

void Eig::arnoldiStep(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
		       Eigen::MatrixXcd &upperHessEigen,
		       field<Complex> *r, double &beta, int j) {

  //%---------------------------------------------------%
  //| STEP 1: Check if the B norm of j-th residual      |
  //| vector is zero. Equivalent to determine whether   |
  //| an exact j-step Arnoldi factorization is present. |
  //%---------------------------------------------------%
  beta = sqrt(blas::norm2(r));
  
  //%--------------------------------%
  //| STEP 2:  v_{j} = r_{j-1}/rnorm |
  //%--------------------------------%  
  blas::ax(1.0/beta, r);
  blas::copy(kSpace[j], r);
  
  //%----------------------------%
  //| STEP 3:  r_{j} = OP*v_{j}; |
  //%----------------------------%  
  OPERATOR(r, kSpace[j], gauge);
  
  //%-------------------------------------%
  //| The following is needed for STEP 5. |
  //| Compute the B-norm of OP*v_{j}.     |
  //%-------------------------------------%

  double wnorm = sqrt(blas::norm2(r));

  //%-----------------------------------------%
  //| Compute the j-th residual corresponding |
  //| to the j step factorization.            |
  //| Use Classical Gram Schmidt and compute: |
  //| w_{j} <-  V_{j}^T * B * OP * v_{j}      |
  //| r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
  //%-----------------------------------------%

  //%------------------------------------------%
  //| Compute the j Fourier coefficients w_{j} |
  //| WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
  //%------------------------------------------%
  //H_{j,i}_j = v_i^dag * r
  for (int i = 0; i < j+1; i++) {
    upperHessEigen(i,j) = blas::cDotProd(kSpace[i], r);
  }
  
  //%--------------------------------------%
  //| Orthogonalize r_{j} against V_{j}.   |
  //| RESID contains OP*v_{j}. See STEP 3. | 
  //%--------------------------------------%
  //r = r - H_{j,i} * v_j 
  for (int i = 0; i < j+1; i++) {
    blas::axpy(-1.0*upperHessEigen(i,j), kSpace[i], r);
  }
  
  if(j > 0) upperHessEigen(j,j-1) = beta;
  
  beta = sqrt(blas::norm2(r));
  
  //%-----------------------------------------------------------%
  //| STEP 5: Re-orthogonalization / Iterative refinement phase |
  //| Maximum NITER_ITREF tries.                                |
  //|                                                           |
  //|          s      = V_{j}^T * B * r_{j}                     |
  //|          r_{j}  = r_{j} - V_{j}*s                         |
  //|          alphaj = alphaj + s_{j}                          |
  //|                                                           |
  //| The stopping criteria used for iterative refinement is    |
  //| discussed in Parlett's book SEP, page 107 and in Gragg &  |
  //| Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
  //| Determine if we need to correct the residual. The goal is |
  //| to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
  //| The following test determines whether the sine of the     |
  //| angle between  OP*x and the computed residual is less     |
  //| than or equal to 0.717.                                   |
  //%-----------------------------------------------------------%

  int orth_iter = 0;
  int orth_iter_max = 100;
  while(beta < 0.717*wnorm && orth_iter < orth_iter_max) {
    
    //%---------------------------------------------------%
    //| Enter the Iterative refinement phase. If further  |
    //| refinement is necessary, loop back here. The loop |
    //| variable is ITER. Perform a step of Classical     |
    //| Gram-Schmidt using all the Arnoldi vectors V_{j}  |
    //%---------------------------------------------------%

    //%---------------------------------------------%
    //| Compute the correction to the residual:     |
    //| r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
    //| The correction to H is v(:,1:J)*H(1:J,1:J)  |
    //| + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
    //%---------------------------------------------%

    wnorm = beta;
    
    // reorthogonalise r against the K space
    //printf("beta = %e < %e: Reorthogonalise at step %d, iter %d\n", beta, 0.717*wnorm, j, orth_iter);
    std::vector<Complex> alpha(j+1, 0.0);
    for(int i=0; i < j+1; i++) {
      alpha[i] = blas::cDotProd(kSpace[i], r);
      upperHessEigen(i,j) += alpha[i];
    }
    
    for(int i=0; i < j+1; i++) {
      blas::axpy(-1.0*alpha[i], kSpace[i], r);
    }
    
    beta = sqrt(blas::norm2(r));
    orth_iter++;
  }
  
  if(orth_iter == orth_iter_max) {
    //%---------------------------------------%
    //| RESID is numerically in the span of V |
    //%---------------------------------------%
    cout << "Eig: Unable to orthonormalise r" << endl;
    exit(0);
  }
}

void Eig::reorder(std::vector<field<Complex> *> kSpace, std::vector<Complex> evals, std::vector<double> residua, int n_kr, Spectrum spectrum) {

  int n = n_kr;
  std::vector<std::tuple<Complex, double, field<Complex>* >> array(n);
  for(int i=0; i<n; i++) array[i] = std::make_tuple(evals[i], residua[i], kSpace[i]);
  
  switch(spectrum) {
  case LM:
    std::sort(array.begin(), array.begin() + n,
	      [] (const std::tuple<Complex, double, field<Complex>*> &a,
		  const std::tuple<Complex, double, field<Complex>*> &b) {
		return (abs(std::get<0>(a)) > abs(std::get<0>(b))); } );
    break;
  case SM:
    std::sort(array.begin(), array.begin() + n,
	      [] (const std::tuple<Complex, double, field<Complex>*> &a,
		  const std::tuple<Complex, double, field<Complex>*> &b) {
		return (abs(std::get<0>(a)) < abs(std::get<0>(b))); } );
    break;
  case LR:
    std::sort(array.begin(), array.begin() + n,
	      [] (const std::tuple<Complex, double, field<Complex>*> &a,
		  const std::tuple<Complex, double, field<Complex>*> &b) {
		return (std::get<0>(a).real() > std::get<0>(b).real()); } );
    break;
  case SR:
    std::sort(array.begin(), array.begin() + n,
	      [] (const std::tuple<Complex, double, field<Complex>*> &a,
		  const std::tuple<Complex, double, field<Complex>*> &b) {
		return (std::get<0>(a).real() < std::get<0>(b).real()); } );
    break;
  case LI:
    std::sort(array.begin(), array.begin() + n,
	      [] (const std::tuple<Complex, double, field<Complex>*> &a,
		  const std::tuple<Complex, double, field<Complex>*> &b) {
		return (std::get<0>(a).imag() > std::get<0>(b).imag()); } );
    break;
  case SI:
    std::sort(array.begin(), array.begin() + n,
	      [] (const std::tuple<Complex, double, field<Complex>*> &a,
		  const std::tuple<Complex, double, field<Complex>*> &b) {
		return (std::get<0>(a).imag() < std::get<0>(b).imag()); } );
    break;
  default: printf("Eig: Undefined spectrum type %d given", spectrum);
    exit(0);
  }
  
  // Repopulate arrays with sorted elements
  for(int i=0; i<n; i++) {
    std::swap(evals[i], std::get<0>(array[i]));
    std::swap(residua[i], std::get<1>(array[i]));
    std::swap(kSpace[i], std::get<2>(array[i]));
  }
}

void Eig::eigensolveFromUpperHess(MatrixXcd &upperHessEigen, MatrixXcd &Qmat,
				   std::vector<Complex> &evals,
				   std::vector<double> &residua,
				   const double beta, int n_kr)
{
  // QR the upper Hessenberg matrix
  Eigen::ComplexSchur<MatrixXcd> schurUH;
  Qmat.setIdentity();
  schurUH.computeFromHessenberg(upperHessEigen, Qmat);
  
  // Extract the upper triangular matrix, eigensolve, then
  // get the eigenvectors of the upper Hessenberg
  MatrixXcd matUpper = MatrixXcd::Zero(n_kr, n_kr);
  matUpper = schurUH.matrixT().triangularView<Eigen::Upper>();
  matUpper.conservativeResize(n_kr, n_kr);
  Eigen::ComplexEigenSolver<MatrixXcd> eigenSolver(matUpper);
  Qmat = schurUH.matrixU() * eigenSolver.eigenvectors();
  
  // Update eigenvalues, residiua, and the Q matrix
  for(int i=0; i<n_kr; i++) {
    evals[i] = eigenSolver.eigenvalues()[i];
    residua[i] = abs(beta * Qmat.col(i)[n_kr-1]);
  }  
}

void Eig::qriteration(MatrixXcd &Rmat, MatrixXcd &Qmat, const int n_kr, const double tol)
{  
  Complex T11, T12, T21, T22, U1, U2;
  double dV;

  //double tol = 1e-15;
  
  // Allocate the rotation matrices.
  std::vector<Complex> R11(n_kr-1);
  std::vector<Complex> R12(n_kr-1);
  std::vector<Complex> R21(n_kr-1);
  std::vector<Complex> R22(n_kr-1);

  // First pass, determine the matrices and do H -> R.
  for(int i = 0; i < n_kr-1; i++) {
    if (abs(Rmat(i+1, i)) < tol) {
      R11[i] = R12[i] = R21[i] = R22[i] = 0.0;
      Rmat(i+1, i) = 0.0;
      continue;
    }
    
    dV = sqrt(norm(Rmat(i, i)) + norm(Rmat(i+1, i)));
    U1 = Rmat(i, i);
    dV = (U1.real() > 0.0) ? dV : -dV;
    U1 += dV;
    U2 = Rmat(i+1, i);
        
    T11 = conj(U1) / dV;
    R11[i] = conj(T11);

    T12 = conj(U2) / dV;
    R12[i] = conj(T12);
    
    T21 = conj(T12) * conj(U1) / U1;
    R21[i] = conj(T21);

    T22 = T12 * U2 / U1;
    R22[i] = conj(T22);

    Rmat(i, i) -= (T11 * Rmat(i, i) + T12 * Rmat(i+1, i));
    Rmat(i+1, i) = 0;
    
#pragma omp parallel for schedule(static,32)  
    for(int j=i+1; j < n_kr; j++) {
      Complex temp = Rmat(i, j);
      Rmat(i, j)   -= (T11 * Rmat(i, j) + T12 * Rmat(i+1, j));      
      Rmat(i+1, j) -= (T21 * temp + T22 * Rmat(i+1, j));
    }
  }

  for(int j = 0; j < n_kr - 1; j++) {
    if(abs(R11[j]) > tol) {
#pragma omp parallel 
      {
#pragma omp for schedule(static,32)  
	for(int i = 0; i < j+2; i++) {
	  Complex temp = Rmat(i, j);
	  Rmat(i, j) -= (R11[j] * temp + R12[j] * Rmat(i, j+1));	
	  Rmat(i, j+1) -= (R21[j] * temp + R22[j] * Rmat(i, j+1));	
	}
	
#pragma omp for schedule(static,32) 
	for(int i = 0; i < n_kr; i++) {
	  Complex temp = Qmat(i, j);
	  Qmat(i, j) -= (R11[j] * temp + R12[j] * Qmat(i, j+1));
	  Qmat(i, j+1) -= (R21[j] * temp + R22[j] * Qmat(i, j+1));	
	}
      }
    }  
  }
}

int Eig::qrFromUpperHess(MatrixXcd &upperHess, MatrixXcd &Qmat, std::vector<Complex> &evals, std::vector<double> &residua, const double beta, const int n_kr, const double tol)
{

  MatrixXcd Rmat = MatrixXcd::Zero(n_kr, n_kr);
  Rmat = upperHess;
  
  //double tol = 1e-15;
  Complex temp, discriminant, sol1, sol2, eval;
  int max_iter = 100000;
  int iter = 0;
  
  for (int i = n_kr-2; i >= 0; i--) {    
    while (iter < max_iter) {
      if(abs(Rmat(i+1, i)) < tol) {
	Rmat(i+1, i) = 0.0;
	break;
      } else {
      
	// Compute the 2 eigenvalues via the quadratic formula
	//----------------------------------------------------

	// The discriminant
	temp = (Rmat(i, i) - Rmat(i+1, i+1)) * (Rmat(i, i) - Rmat(i+1, i+1)) / 4.0;
	discriminant = sqrt(Rmat(i+1, i) * Rmat(i, i+1) + temp);

	// Reuse temp
	temp = (Rmat(i, i) + Rmat(i+1, i+1))/2.0;
	
	sol1 = temp - Rmat(i+1, i+1) + discriminant;
	sol2 = temp - Rmat(i+1, i+1) - discriminant;

	// Deduce the better eval to shift
	eval = Rmat(i+1, i+1) + (norm(sol1) < norm(sol2) ? sol1 : sol2);
	
	// Shift the eigenvalue
	for(int j = 0; j < n_kr; j++) Rmat(j, j) -= eval;
	
	// Do the QR iteration
	qriteration(Rmat, Qmat, n_kr, tol);
	
	// Shift back
	for(int j = 0; j < n_kr; j++) Rmat(j, j) += eval;
      }
      iter++;    
    } 
  }

  MatrixXcd matUpper = MatrixXcd::Zero(n_kr, n_kr);
  matUpper = Rmat.triangularView<Eigen::Upper>();
  matUpper.conservativeResize(n_kr, n_kr);
  Eigen::ComplexEigenSolver<MatrixXcd> eigenSolver(matUpper);
  Qmat *= eigenSolver.eigenvectors();
  
  // Update eigenvalues, residiua, and the Q matrix
  for(int i=0; i<n_kr; i++) {
    evals[i] = eigenSolver.eigenvalues()[i];
    residua[i] = abs(beta * Qmat.col(i)[n_kr-1]);
  }
  
  if(eig_verbose) printf("Eig: eigensystem iterations = %d\n", iter);

  return iter;  
}

void Eig::rotateVecsComplex(std::vector<field<Complex> *> vecs, Eigen::MatrixXcd mat, int num_locked, int iter_keep, int n_kr) {

  //loop over rows of V_k
  int Nx = vecs[0]->p.Nx;
  int Ny = vecs[0]->p.Ny;
#pragma omp parallel for
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++) {
	
	//put jth row of V_k in temp
	Complex tmp[n_kr];  
	for(int i=0; i<n_kr; i++) {
	  tmp[i] = vecs[i+num_locked]->read(x,y,mu);   
	}
	
	//take product of jth row of V_k and ith column of mat (ith eigenvector of T_k)
	for(int i=0; i<iter_keep; i++) {
	  Complex sum = 0.0;
	  //Loop over elements to get the y_i[j]th element 
	  for(int l=0; l<n_kr; l++) {
	    sum += tmp[l]*mat.col(i)[l];
	  }
	  
	  //Update the Ritz vector
	  vecs[i+num_locked]->write(x,y,mu, sum);
	}
      }
}

void Eig::computeEvals(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace, std::vector<double> &residua, std::vector<Complex> &evals, int n_ev) {
  
  //temp vector
  field<Complex> *temp = new field<Complex>(gauge->p);
  for (int i = 0; i < n_ev; i++) {
    
    // r = A * v_i
    OPERATOR(temp, kSpace[i], gauge);

    // lambda_i = v_i^dag A v_i / (v_i^dag * v_i)
    evals[i] = blas::cDotProd(kSpace[i], temp);

    // Measure ||lambda_i*v_i - A*v_i||
    Complex n_unit(-1.0, 0.0);
    blas::axpby(evals[i], kSpace[i], n_unit, temp);
    residua[i] = sqrt(blas::norm2(temp));
  }
  for (int i = 0; i < n_conv; i++) {
    double norm = sqrt(blas::norm2(kSpace[i]));
    if (residua[i] > 10*tol) printf("Eig CHECK FAIL: EigValue[%04d]:||(%+.8e, %+.8e)|| = %+.8e residual %.8e norm = %.8e\n", i, evals[i].real(), evals[i].imag(), abs(evals[i]), residua[i], norm);
  }
  delete temp;
}

void Eig::solve(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
		std::vector<Complex> &evals) {
  if(gauge->p.use_feast) feast(gauge, kSpace, evals);
  else                   iram(gauge, kSpace, evals);  
}

void Eig::feast(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
		std::vector<Complex> &evals) {

  struct timeval start, end, total_start, total_end;
  gettimeofday(&total_start, NULL);  
  // timing variables
  double t_total =  0.0;
  double t_init = 0;
  double t_crs = 0.0;
  double t_feast = 0.0;

  bool converged = false;
  
#ifdef USE_FEAST
  // START init
  //---------------------------------------------------------  
  gettimeofday(&start, NULL);  
    
  // Compute CRS matrix for FEAST
  // Fortran indexing applies... be careful!
  std::vector<double> crs_elems;
  std::vector<int> crs_col_idx;
  std::vector<int> crs_row_idx;

  string name;
  name = "data/eig/crs_LX" + to_string(Nx) + "_LY" + to_string(Ny) + ".dat";
  fstream input_file;
  input_file.open(name);
  string val;
  if(input_file.is_open()) {
    if(verbosity) cout << "IO: CRS wisdom found for " << name << endl;
    getline(input_file, val);
    int n_row = 2*Nx*Ny + 1;
    int n_col = stoi(val);
    crs_col_idx.reserve(n_col);
    for(int i=0; i<n_col; i++) {
      getline(input_file, val);
      crs_col_idx.push_back(stoi(val));    
    }
      
    for(int i=0; i<n_row; i++) {
      getline(input_file, val);
      crs_row_idx.push_back(stoi(val));    
    }

    crs_elems.reserve(crs_row_idx[n_row-1]);
    // Get the ith column data, conjugate it, it becomes the
    // the data on the ith row. 
    field<Complex> *op_column = new field<Complex>(gauge->p);
    field<Complex> *point_source = new field<Complex>(gauge->p);
    int col_idx = 0;
    for(int i=0; i<2*Nx*Ny; i++) {
      blas::zero(point_source);
      point_source->write(i, cUnit);
      OPERATOR(op_column, point_source, gauge);
      
      // Store the current row index
      // Loop over this row, get non-zero data
      for(int j=0; j<crs_row_idx[i+1] - crs_row_idx[i]; j++) {      	  
	// Store data in the array
	crs_elems.push_back( op_column->elem(crs_col_idx[col_idx] - 1).real());
	crs_elems.push_back(-op_column->elem(crs_col_idx[col_idx] - 1).imag());
	col_idx++;
      }
    }

    delete op_column;
    delete point_source;

  } else {
    if(verbosity) cout << "IO: NO CRS wisdom found for " << name << endl;
    // Row index always starts at 1
    int current_crs_row_idx = 1;  

    field<Complex> *op_column = new field<Complex>(gauge->p);
    field<Complex> *point_source = new field<Complex>(gauge->p);
      
    for(int i=0; i<2*Nx*Ny; i++) {      
      // Get the ith column data, conjugate it, it becomes the
      // the data on the ith row. 
      blas::zero(point_source);
      point_source->write(i, cUnit);
      OPERATOR(op_column, point_source, gauge);
      
      // Store the current row index
      crs_row_idx.push_back(current_crs_row_idx);
      // Loop over this row, get non-zero data
      for(int j=0; j<2*Nx*Ny; j++) {      
	if(fabs(op_column->elem(j).real()) >= 1e-12 || fabs(op_column->elem(j).imag()) > 1e-12) {
	  
	  // Store data in the array
	  crs_elems.push_back( op_column->elem(j).real());
	  crs_elems.push_back(-op_column->elem(j).imag());
	  
	  // Store the column index of this row
	  crs_col_idx.push_back(j+1);
	  
	  // Increment the row index
	  current_crs_row_idx++;
	}
      }
    }
    delete op_column;
    delete point_source;
    crs_row_idx.push_back(current_crs_row_idx);
    
    cout << "IO: Writing CRS data " << name << endl;
    fstream output_file;
    output_file.open(name,ios::in|ios::out|ios::trunc);  
    
    output_file << crs_col_idx.size() << endl;
    for(int i=0; i<crs_col_idx.size(); i++) output_file << crs_col_idx[i] << endl;
    for(int i=0; i<crs_row_idx.size(); i++) output_file << crs_row_idx[i] << endl;
  }
    
  gettimeofday(&end, NULL);  
  t_crs += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
  if(verbosity) {
    cout << "FEAST: CRS construction time = " << t_crs << endl;
    cout << "FEAST: CRS elems = " << crs_elems.size()/2 << endl;
    cout << "FEAST: CRS row idx size = " << crs_row_idx.size() << " Should equal " << 2*Nx*Ny + 1<< endl;
    //for(int g=0; g<crs_row_idx.size(); g++) cout << "crs_row_idx["<<g<<"] = " << crs_row_idx[g] << endl;
    //for(int g=0; g<crs_col_idx.size(); g++) cout << "crs_col_idx["<<g<<"] = " << crs_col_idx[g] << " elem: " << crs_elems[2*g] << " " << crs_elems[2*g+1] << endl;
  }

  // START init
  gettimeofday(&start, NULL);  

  int  M = n_conv;
  char UPLO='F';
    
  int N  = 2*Nx*Ny;
  int M0 = feast_M0; // M0 >= M
    
  int fpm[64];
  double epsout;
  int loop; 
  int info;
  double Emid[2];
  Emid[0] = 0.0;
  Emid[1] = feast_Emax;;

  // Extend FEAST space and evals
  if(inspection_counter == 0) {
    kSpace_feast.resize(M0);
    for(int i=0; i<M0; i++) kSpace_feast[i] = new field<Complex>(gauge->p);
    evals_feast.resize(M0);
  }
    
  // FEAST defaults
  feastinit(fpm);  

  // Adjust FEAST params for this problem
  fpm[0] = eig_verbose ? 1 : 0;  // FEAST verbosity
  fpm[1] = feast_Ncontour;                  // FEAST contour points
  fpm[2] = (int)(-log10(tol)); // FEAST tol
  fpm[3] = max_restarts;
  fpm[4] = inspection_counter > 0 ? 1 : 0;
  fpm[9] = 0;
  fpm[13] = 0;
  fpm[15] = 0; // gauss/trap/zolatorov
  fpm[39] = 0; // Provide search interval Emid[0] -> Emid[1]
  fpm[42] = 0; // IFEAST
    
  std::vector<double> E(M0, 0.0);     // eigenvalues
  std::vector<double> res(M0, 0.0);   // residua (if needed)
  std::vector<double> X(2*N*M0, 0.0); // eigenvectors:factor 2 for complex
  if(fpm[4] == 1) {
    // Transfer data from krylov space to FEAST
    Complex elem;      
    for(int vec=0; vec<M0; vec++) {
      E[vec] = evals_feast[vec].real();
      for(int x=0; x<Nx; x++) {
	for(int y=0; y<Ny; y++) {
	  for(int mu=0; mu<2; mu++) {
	    int glo_idx = 2*((y * gauge->p.Nx) + x) + mu; // Global index
	      
	    elem = kSpace_feast[vec]->read(x, y, mu);
	    X[2*(glo_idx + N*vec)    ] = elem.real();
	    X[2*(glo_idx + N*vec) + 1] = elem.imag();
	  }
	}
      }
    }
  }
  gettimeofday(&end, NULL);  
  t_init += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

  // FEST SOLVER
  gettimeofday(&start, NULL);  
  zfeast_hcsrev(&UPLO, &N, crs_elems.data(), crs_row_idx.data(), crs_col_idx.data(), fpm, &epsout, &loop, &Emid[0], &Emid[1], &M0, E.data(), X.data(), &M, res.data(), &info);
  gettimeofday(&end, NULL);  
  t_feast += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
  if (info!=0) {
    // Post computation report  
    printf("FEAST ERROR: Failed to compute the requested %d vectors with with a %d search space and a %d Krylov space: error %d\n", n_conv, M, M0, info);
    exit(0);
  } else {

    converged = true;
    
    // FINISH init      
    gettimeofday(&start, NULL);
      
    // Migrate data to kSpace space
    Complex elem;
      
    for(int vec=0; vec<n_conv; vec++) {
      evals[vec].real(E[vec]);
      evals[vec].imag(0.0);
      if(deflationSpaceExists()) {
	evals_defl[vec].real(E[vec]);
	evals_defl[vec].imag(0.0);
      }
	
      for(int x=0; x<Nx; x++) {
	for(int y=0; y<Ny; y++) {
	  for(int mu=0; mu<2; mu++) {
	    int glo_idx = 2*((y * gauge->p.Nx) + x) + mu; // Global index
	      
	    elem.real(X[2*(glo_idx   + N*vec)]);
	    elem.imag(X[2*(glo_idx   + N*vec)+1]);	
	    kSpace[vec]->write(glo_idx, elem);
	    if(deflationSpaceExists()) kSpace_defl[vec]->write(glo_idx, elem);
	  }
	}
      }
    }

    if(kSpace_feast.size() > 0) {
      // Migrate data to FEAST space
      for(int vec=0; vec<M0; vec++) {
	evals_feast[vec].real(E[vec]);
	evals_feast[vec].imag(0.0);
	for(int x=0; x<Nx; x++) {
	  for(int y=0; y<Ny; y++) {
	    for(int mu=0; mu<2; mu++) {
	      int glo_idx = 2*((y * gauge->p.Nx) + x) + mu; // Global index
		
	      elem.real(X[2*(glo_idx   + N*vec)  ]);
	      elem.imag(X[2*(glo_idx   + N*vec)+1]);	
	      kSpace_feast[vec]->write(glo_idx, elem);
	    }
	  }
	}
      }
    }

    if(eig_verbose) {
      // Measure the kSpace inner products
      std::vector<Complex> inner_products(n_conv * n_conv, 0.0);
      for(int i=0; i<n_conv; i++) {
	for(int j=0; j<n_conv; j++) {
	  inner_products[i + n_conv * j] = blas::cDotProd(kSpace[j], kSpace[i]);
	  if(i==j && 1 - abs(inner_products[i + n_conv * j]) > 1e-12) cout << "FEAST ortho fail at (" <<i<< "," <<j<< ")" << endl;
	  if(i!=j && abs(inner_products[i + n_conv * j]) > 1e-12) cout << "FEAST ortho fail at (" <<i<< "," <<j<< ")" << endl;
	}
      }
    }

    gettimeofday(&end, NULL);  
    t_init += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

    gettimeofday(&total_end, NULL);  
    t_total += ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
      
    if(verbosity) {
      printf("FEAST: Computed the requested %d vectors with a %d search space and a %d Krylov space\n", n_conv, M, M0);    
      for (int i = 0; i < n_conv; i++) {
	printf("FEAST: EigValue[%04d]:||(%+.8e, %+.8e)|| = %+.8e residual %.8e\n", i, evals[i].real(), evals[i].imag(), abs(evals[i]), res[i]);
      }
	
      cout << "FEAST: Timings:" << endl;
      cout << "FEAST: crs = " << t_crs << endl;
      cout << "FEAST: init = " << t_init << endl;
      cout << "FEAST: compute = " << t_feast << endl;
      double t_sum = t_feast + t_init + t_crs;
      cout << "FEAST: missing = " << t_total << " - " << t_sum << " = " << t_total - t_sum << " ("<<100*(t_total - t_sum)/t_total<<"%)" << endl;  
    }      
  }
  
  // Clean up
  kSpace.resize(n_conv);
  evals.resize(n_conv);
  
#else
  cout << "ERROR: FEAST not installed " << endl;
  exit(0);
#endif
}
  
void Eig::iram(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
	       std::vector<Complex> &evals) {
  
  struct timeval start, end, total_start, total_end;
  gettimeofday(&total_start, NULL);
  
  // timing variables
  double t_total =  0.0;
  double t_init = 0;
  double t_sort = 0;
  double t_eigen = 0.0;
  double t_compute = 0.0;
  double t_QR = 0.0;
  double t_EV = 0.0;  

  bool converged = false;
    
  if (!(n_kr > n_ev + 6) && n_kr != (gauge->p.Nx * gauge->p.Ny * 2)) {
    printf("IRAM WARNING: n_kr=%d must be greater than n_ev+6=%d\n", n_kr, n_ev + 6);
    exit(0);
  }

  if(eig_verbose) {
    printf("IRAM: n_kr = %d\n", n_kr);
    printf("IRAM: n_ev = %d\n", n_ev);
    printf("IRAM: n_conv = %d\n", n_conv);
    printf("IRAM: Restarts = %d\n", max_restarts);
    printf("IRAM: tol = %e\n", tol);
  }

  // Construct objects for IRAM.
  //---------------------------------------------------------------------
  //Eigenvalues and their residua
  std::vector<double> residua(n_kr, 0.0);
    
  // Upper Hessenberg matrix
  MatrixXcd upperHessEigen = MatrixXcd::Zero(n_kr, n_kr);
    
  // Residual vector. Also used as temp vector(s)
  field<Complex> *r = new field<Complex>(gauge->p);
    
  // Extend Krylov space and evals
  kSpace.resize(n_kr);
  for(int i=n_conv; i<n_kr; i++) kSpace[i] = new field<Complex>(gauge->p);
  evals.resize(n_kr);
    
  // Eigen objects for Arnoldi vector rotation and QR shifts
  MatrixXcd Qmat = MatrixXcd::Identity(n_kr, n_kr);
  MatrixXcd sigma = MatrixXcd::Identity(n_kr, n_kr);
    
  double epsilon = 2e-16;
  double epsilon23 = pow(epsilon, 2.0/3.0);
  double beta = 0.0;
  int iter = 0;
  int restart_iter = 0;
  int iter_converged = 0;
  int iter_keep = 0;
  int num_converged = 0;
  int num_keep = 0;
  // END init
  //---------------------------------------------------------
  gettimeofday(&end, NULL);  
  t_init += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
  // START compute 
  //---------------------------------------------------------
  gettimeofday(&start, NULL);  
  // Populate source with randoms.
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	Complex temp(drand48(), drand48());
	r->write(x,y,mu, temp);  
      }
    }
  }
    
  //Place initial source in range of mat
  OPERATOR(kSpace[0], r, gauge);
  blas::copy(r, kSpace[0]);
    
  // START Eig
  // Implicitly restarted Arnoldi method for asymmetric eigenvalue problems
  if(eig_verbose) printf("IRAM: Start IRAM solution\n");
  //----------------------------------------------------------------------
    
  // Loop over restart iterations.
  num_keep = 0;
  while(restart_iter < max_restarts && !converged) {
    gettimeofday(&start, NULL); 
    for (int step = num_keep; step < n_kr; step++) arnoldiStep(gauge, kSpace, upperHessEigen, r, beta, step);
    iter += (n_kr - num_keep);
    gettimeofday(&end, NULL);  
    t_compute += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
      
    // Compute Ritz and bounds
    gettimeofday(&start, NULL);
    qrFromUpperHess(upperHessEigen, Qmat, evals, residua, beta, n_kr, tol/10);
    gettimeofday(&end, NULL);  
    t_EV += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
      
    num_keep = n_ev;    
    int nshifts = n_kr - num_keep;
      
    // Emulate zngets: sort the unwanted Ritz to the start of the arrays, then
    // sort the first (n_kr - n_ev) bounds to be first for forward stability
    gettimeofday(&start, NULL); 
    // Sort to put unwanted Ritz(evals) first
    zsortc(spectrum, n_kr, evals, residua);    
    // Sort to put smallest Ritz errors(residua) first
    zsortc(LM, nshifts, residua, evals);
      
    // Convergence test
    iter_converged = 0;
    for(int i=0; i<n_ev; i++) {
      int idx = n_kr - 1 - i;
      double rtemp = std::max(epsilon23, abs(evals[idx]));
      if(residua[idx] < tol * rtemp) {
	iter_converged++;
	if(eig_verbose) printf("IRAM: residuum[%d] = %e, cond = %e\n", i, residua[idx], tol * abs(evals[idx]));
      } else {
	break;
      }
    }    
      
    //%---------------------------------------------------------%
    //| Count the number of unwanted Ritz values that have zero |
    //| Ritz estimates. If any Ritz estimates are equal to zero |
    //| then a leading block of H of order equal to at least    |
    //| the number of Ritz values with zero Ritz estimates has  |
    //| split off. None of these Ritz values may be removed by  |
    //| shifting. Decrease NP the number of shifts to apply. If |
    //| no shifts may be applied, then prepare to exit          |
    //%---------------------------------------------------------%
      
    int num_keep0 = num_keep;
    iter_keep = std::min(iter_converged + (n_kr - num_converged) / 2, n_kr - 12);
      
    num_converged = iter_converged;
    num_keep = iter_keep;      
    nshifts = n_kr - num_keep;
      
    if(verbosity) printf("IRAM: %04d converged eigenvalues at iter %d\n", num_converged, restart_iter);
      
    int nshifts0 = nshifts;
    for(int i=0; i<nshifts0; i++) {
      if(residua[i] <= epsilon) {
	nshifts--;
      }
    }
      
    if(nshifts == 0 && num_converged < n_ev) {
      cout << "IRAM ERROR: No shifts can be applied" << endl;
      exit(0);
    }    
      
    gettimeofday(&end, NULL);
    t_sort += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
      
    if (num_converged >= n_conv) {
      converged = true;
	
      // Compute Eigenvalues
      gettimeofday(&start, NULL);
      Qmat.setIdentity();
      qrFromUpperHess(upperHessEigen, Qmat, evals, residua, beta, n_kr, 1e-15);
      gettimeofday(&end, NULL);  
      t_EV += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	
      gettimeofday(&start, NULL);
      rotateVecsComplex(kSpace, Qmat, 0, n_kr, n_kr);
      reorder(kSpace, evals, residua, n_kr, spectrum);
      computeEvals(gauge, kSpace, residua, evals, n_kr);
      gettimeofday(&end, NULL);  
      t_compute += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	
    } else if (restart_iter < max_restarts) {
	
      //%-------------------------------------------------%
      //| Do not have all the requested eigenvalues yet.  |
      //| To prevent possible stagnation, adjust the size |
      //| of NEV.                                         |
      //| If the size of NEV was just increased resort    |
      //| the eigenvalues.                                |
      //%-------------------------------------------------%
	
      if(num_keep0 < num_keep) {
	gettimeofday(&start, NULL); 
	// Emulate zngets: sort the unwanted Ritz to the start of the arrays, then
	// sort the first (n_kr - n_ev) bounds to be first for forward stability	
	// Sort to put unwanted Ritz(evals) first
	zsortc(spectrum, n_kr, evals, residua);
	// Sort to put smallest Ritz errors(residua) first
	zsortc(LM, nshifts, residua, evals);
	gettimeofday(&end, NULL);  
	t_sort += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
      }
	
      //%---------------------------------------------------------%
      //| Apply the NP implicit shifts by QR bulge chasing.       |
      //| Each shift is applied to the whole upper Hessenberg     |
      //| matrix H.                                               |
      //%---------------------------------------------------------%
	
      gettimeofday(&start, NULL); 
      Qmat.setIdentity();
      sigma.setIdentity();      
      for(int i=0; i<nshifts; i++){	
	sigma.setIdentity();
	sigma *= evals[i];
	upperHessEigen -= sigma;
	qriteration(upperHessEigen, Qmat, n_kr, tol/10);	
	upperHessEigen += sigma;	
      }
      gettimeofday(&end, NULL);  
      t_QR += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	
      gettimeofday(&start, NULL); 
      rotateVecsComplex(kSpace, Qmat, 0, num_keep+1, n_kr);
	
      //%-------------------------------------%
      //| Update the residual vector:         |
      //|    r <- sigmak*r + betak*v(:,kev+1) |
      //| where                               |
      //|    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
      //|    betak = e_{kev+1}'*H*e_{kev}     |
      //%-------------------------------------%
	
      blas::axpby(upperHessEigen(num_keep, num_keep-1), kSpace[num_keep], Qmat(n_kr-1, num_keep-1), r);
      gettimeofday(&end, NULL);  
      t_compute += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	
      if(sqrt(blas::norm2(r)) < epsilon) {
	printf("IRAM ERROR: Congratulations! You have reached an invariant subspace at iter %d, beta = %e\n", restart_iter, sqrt(blas::norm2(r)));
	exit(0);
      }
    }
    restart_iter++;    
  }
    
  gettimeofday(&total_end, NULL);  
  t_total = ((total_end.tv_sec - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
  
  // Post computation report  
  if (!converged) {    
    printf("IRAM ERROR: Failed to compute the requested %d vectors with with a %d search space and a %d Krylov space in %d restart_steps and %d OPs.\n", n_conv, n_ev, n_kr, restart_iter, iter);
    exit(0);
  } else {
    if(verbosity) {
      printf("IRAM: Computed the requested %d vectors with a %d search space and a %d Krylov space in %d restart_steps and %d OPs in %e secs.\n", n_conv, n_ev, n_kr, restart_iter, iter, (t_compute + t_sort + t_EV + t_QR));    
      for (int i = 0; i < n_conv; i++) {
	printf("IRAM: EigValue[%04d]:||(%+.8e, %+.8e)|| = %+.8e residual %.8e\n", i, evals[i].real(), evals[i].imag(), abs(evals[i]), residua[i]);
      }
    }    
  }  

  if(eig_verbose) {
    // Measure the kSpace inner products
    std::vector<Complex> inner_products(n_conv * n_conv, 0.0);
    for(int i=0; i<n_conv; i++) {
      for(int j=0; j<n_conv; j++) {
	inner_products[i + n_conv * j] = blas::cDotProd(kSpace[j], kSpace[i]);
	if(i==j && (1 - abs(inner_products[i + n_conv * j]) > 1e-12)) cout << "IRAM ortho fail at (" <<i<< "," <<j<< "): " << 1 - abs(inner_products[i + n_conv * j]) << endl;
	if(i!=j && (abs(inner_products[i + n_conv * j]) > 1e-12)) cout << "IRAM ortho fail at (" <<i<< "," <<j<< "): " << abs(inner_products[i + n_conv * j]) << endl;
      }
    }
  }
  if(verbosity) {
    cout << "IRAM: Timings:" << endl;
    cout << "IRAM: init = " << t_init << endl;
    cout << "IRAM: compute = " << t_compute << endl;
    cout << "IRAM: sort = " << t_sort << endl;
    cout << "IRAM: EV = " << t_EV << endl;
    cout << "IRAM: QR = " << t_QR << endl;
    cout << "IRAM: missing = " << (t_total) << " - " << (t_compute + t_init + t_sort + t_EV + t_QR + t_eigen) << " = " << (t_total - (t_compute + t_init + t_sort + t_EV + t_QR + t_eigen)) << " ("<<(100*((t_total - (t_compute + t_init + t_sort + t_EV + t_QR + t_eigen))))/t_total<<"%)" << endl;  
  }
  
  // Clean up
  kSpace.resize(n_conv);
  evals.resize(n_conv);    
}

void Eig::inspectEvolvedSpectrum(const field<Complex> *gauge, int iter) {
  
  std::vector<field<Complex>*> kSpace;
  std::vector<Complex> evals;
  
  std::vector<field<Complex>*> kSpace_aux;
  std::vector<Complex> evals_aux;

  // Allocate space to hold the previous kSpace
  if(kSpace_pre.size() == 0) prepareKrylovSpace(kSpace_pre, evals_pre, gauge->p);  
  
  // Allocate space for this inspection
  prepareKrylovSpace(kSpace, evals, gauge->p);

  // Compute eigenspectrum of this gauge.
  if(inspection_counter == 0) {
    // A deflation space was just constructed, copy it to the
    // kSpace container and save.
    if(eig_verbose) cout << "Eig: First inspection call, copy deflation space" << endl;
    for(int i=0; i<n_conv; i++) {
      blas::copy(kSpace[i], kSpace_defl[i]);
      evals[i] = evals_defl[i];
    }
    
  } else {
    // If inspection_counter > 0, the gauge field has changed and
    // we need to recompute the space for inspecton
    if(eig_verbose) cout << "Eig: Computing new spectrum" << endl;
    solve(gauge, kSpace, evals);
  }
    
  string name;
  char fname[256];
  FILE *fp;
  
  // Dump eigenvalue data
  name = "data/eig/eigenvalues" + to_string(iter);
  constructName(name, gauge->p);
  name += ".dat";
  snprintf(fname, 100, "%s", name.c_str());	  
  fp = fopen(fname, "a");
  for(int i=0; i<n_conv; i++) {
    fprintf(fp, "%d %d %.16e %.16e\n",
	    inspection_counter,
	    i,
	    evals[i].real(),
	    evals[i].imag());
  }
  fclose(fp);

  // Dump eigenvector data
  name = "data/eig/eigenvectors" + to_string(iter);
  constructName(name, gauge->p);
  name += ".dat";
  snprintf(fname, 100, "%s", name.c_str());	  
  fp = fopen(fname, "a");  
  for(int i=0; i<n_conv; i++) {
    // Data to the ith eigenvector
    for(int x=0; x<gauge->p.Nx; x++) {
      for(int y=0; y<gauge->p.Ny; y++) {
	for(int mu=0; mu<2; mu++) {
	  int glo_idx = 2*((y * gauge->p.Nx) + x) + mu; // Global index
	  
	  fprintf(fp, "%d   %d %d %d %d %.16e %.16e %.16e %.16e\n",
		  inspection_counter, i, x, y, mu,
		  kSpace[i]->data[glo_idx].real(),
		  kSpace[i]->data[glo_idx].imag(),
		  abs(kSpace[i]->data[glo_idx]),
		  atan2(kSpace[i]->data[glo_idx].imag(), kSpace[i]->data[glo_idx].real()));
	}
      }
    }
  }
  fclose(fp);
  
  if(inspection_counter > 0) {
    // Measure the kSpace inner products
    std::vector<Complex> inner_products(n_conv * n_conv);
    std::vector<std::vector<int>> swap_mat(n_conv, std::vector<int> (n_conv, 0));
    for(int i=0; i<n_conv; i++) {
      for(int j=0; j<n_conv; j++) {
	inner_products[i + n_conv * j] = blas::cDotProd(kSpace_pre[j], kSpace[i]);
	swap_mat[i][j] = int_round(abs(inner_products[i + n_conv * j]));
      }
    }

#if 0
    for(int i=0; i<n_conv; i++) {
      for(int j=0; j<n_conv; j++) {
	swap_mat[i][j] == 1 ? cout << swap_mat[i][j] << " " : cout << "  ";
      }
      cout << endl;
    }
#endif
    
    // Measure kSpace overlaps
    std::vector<double> overlaps(n_conv, 0.0);
    for(int i=0; i<n_conv; i++) 
      for(int j=0; j<n_conv; j++) overlaps[i] += abs(inner_products[i + n_conv * j] * conj(inner_products[i + n_conv * j]));
    
    // Dump inner product data
    name = "data/eig/inner_products" + to_string(iter);
    constructName(name, gauge->p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());	  
    fp = fopen(fname, "a");
    for(int i=0; i<n_conv; i++) {
      for(int j=0; j<n_conv; j++) {
	fprintf(fp, "%d   %d %d %.16e %.16e %.16e\n",
		inspection_counter, i, j,
		abs(inner_products[i + n_conv * j]),
		inner_products[i + n_conv * j].real(),
		inner_products[i + n_conv * j].imag());
      }
    }
    fprintf(fp,"\n");
    fclose(fp);
    
    // Dump overlap data
    name = "data/eig/overlap" + to_string(iter);
    constructName(name, gauge->p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());	  
    fp = fopen(fname, "a");
    for(int i=0; i<n_conv; i++) {
      fprintf(fp, "%d   %d %.16e\n",
	      inspection_counter, i,
	      overlaps[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
  }

  if(use_compressed_space) {
    // Compute an new MG deflation space 
    if(inspection_counter > 0) {
      prepareKrylovSpace(kSpace_aux, evals_aux, gauge->p);
      computeMGDeflationSpace(kSpace_aux, kSpace, gauge);
    }
    
    // Dump block data
    name = "data/eig/block_data_orth" + to_string(iter);
    constructName(name, gauge->p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());	  
    fp = fopen(fname, "a");
    
    for(int by=0; by<block_scheme[1]; by++) {
      for(int bx=0; bx<block_scheme[0]; bx++) {
	int blk_idx = by * block_scheme[0] + bx;
	for(int i=0; i<n_low; i++) {
	  for(int j=0; j<block_size; j++) {
	    fprintf(fp, "%d   %d %d %d %d %.16e %.16e %.16e %.16e\n",
		    inspection_counter, bx, by, i, j,
		    block_data_ortho[blk_idx][block_size * i + j].real(),
		    block_data_ortho[blk_idx][block_size * i + j].imag(),
		    abs(block_data_ortho[blk_idx][block_size * i + j]),
		    atan2(block_data_ortho[blk_idx][block_size * i + j].imag(), block_data_ortho[blk_idx][block_size * i + j].real()));
	  }
	}    
      }
    }
    fclose(fp);
    
    // Dump block coeffs
    name = "data/eig/block_data_coef" + to_string(iter);
    constructName(name, gauge->p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());	  
    fp = fopen(fname, "a");
    
    for(int by=0; by<block_scheme[1]; by++) {
      for(int bx=0; bx<block_scheme[0]; bx++) {
	int blk_idx = by * block_scheme[0] + bx;
	for(int i=0; i<n_conv; i++) {
	  for(int j=0; j<n_low; j++) {
	    fprintf(fp, "%d   %d %d %d %d %.16e %.16e %.16e %.16e\n",
		    inspection_counter, bx, by, i, j,
		    block_coeffs[blk_idx][n_low * i + j].real(),
		    block_coeffs[blk_idx][n_low * i + j].imag(),
		    abs(block_coeffs[blk_idx][n_low * i + j]),
		    atan2(block_coeffs[blk_idx][n_low * i + j].imag(), block_coeffs[blk_idx][n_low * i + j].real()));
	  }
	}    
      }
    }
    fclose(fp);
  }
  
  // Copy the new space into the old space for overlap comparision.
  for(int i=0; i<n_conv; i++) {
    blas::copy(kSpace_pre[i], kSpace[i]);
    evals_pre[i] = evals[i];
  }
  
  // Clean up
  for(int i=0; i<kSpace.size(); i++) delete kSpace[i];
  for(int i=0; i<kSpace_aux.size(); i++) delete kSpace_aux[i];
  
  inspection_counter++;
}

void Eig::computeDeflationSpace(const field<Complex> *gauge){
  prepareKrylovSpace(kSpace_defl, evals_defl, gauge->p);
  solve(gauge, kSpace_defl, evals_defl);  

  if(use_compressed_space) {
    prepareKrylovSpace(kSpace_mg, evals_mg, gauge->p);
    computeMGDeflationSpace(kSpace_mg, kSpace_defl, gauge);
  }
}

void Eig::computeMGDeflationSpace(std::vector<field<Complex>*> &kSpace_out, const std::vector<field<Complex>*> &kSpace_in,
				   const field<Complex> *gauge){

  struct timeval start, end;
  
  gettimeofday(&start, NULL);    
  // Compress kSpace into block_data_ortho and block_coeffs...
  blockCompress(kSpace_in);
  // ...then expand to into kSpace
  blockExpand(kSpace_out);
  gettimeofday(&end, NULL);  
  
  // Compute the eigenvalues and residua using the reconstructed kSpace
  std::vector<double> resid(n_conv, 0.0);
  computeEvals(gauge, kSpace_out, resid, evals_mg, n_conv);
  
  if(verbosity) {
    for (int i = 0; i < n_conv; i++) {
      printf("Compression: EigValue[%04d]: ||(%+.8e, %+.8e)|| = %+.8e residual %.8e\n", i, evals_mg[i].real(), evals_mg[i].imag(), abs(evals_mg[i]), resid[i]);
    }
    
    // Check compression ratio
    int pre = 2 * n_conv * Nx * Ny;
    int post = n_blocks * n_low * (block_size + n_conv);
    cout << "Compression: Algorithmic compression report: " << endl;
    cout << "Compression: Complex(double) elems pre = " << pre << " Complex(double) elems post = " << post << endl;
    cout << "Compression: Ratio: " << (100.0 * post)/pre << "% of original data." << endl;
    cout << "Compression: " << n_low << " low eigenvectors used " << endl;
    cout << "Compression: " << (n_conv - n_low) << " high eigenvectors reconstructed " << endl;
    cout << "Compression: Compress/decompress time = " << ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 << endl;
  }
}

void Eig::prepareKrylovSpace(std::vector<field<Complex>*> &kSpace,
			      std::vector<Complex> &evals,
			      const Param p) {  
  
  kSpace.reserve(eig_param.n_conv);
  for(int i=0; i<eig_param.n_conv; i++) kSpace.push_back(new field<Complex>(p));
  evals.resize(eig_param.n_conv);
}

// Read the block data from the (iEig)th vector in kSpace 
void Eig::readVectorToBlock(const std::vector<field<Complex>*> &kSpace, std::vector<std::vector<Complex>> &block_data) {
  
  for(int i=0; i<n_conv; i++) {
    for(int by=0; by<block_scheme[1]; by++) {
      for(int bx=0; bx<block_scheme[0]; bx++) {
	int blk_idx = by * block_scheme[0] + bx;
	
	// Location of the start of the desired block (x runs fastest)
	int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));
	
	for(int nx=0; nx<x_block_size; nx++) {
	  for(int ny=0; ny<y_block_size; ny++) {
	    for(int mu=0; mu<2; mu++) {
	      
	      int loc_idx = 2*(nx + x_block_size * ny) + mu;      // Local index
	      int glo_idx = blk_offset + 2*((ny * Nx) + nx) + mu; // Global index
	      block_data[blk_idx][block_size * i + loc_idx] = kSpace[i]->data[glo_idx];
	    }
	  }
	}
      }
    }
  }
}

void Eig::blockCompress(const std::vector<field<Complex>*> &kSpace) {
  
  // Object to hold the blocked eigenvector data
  std::vector<std::vector<Complex>> block_data(n_blocks, std::vector<Complex> (n_conv * block_size, 0.0));
  
  // Copy data from the eigenvector array into a block array
  readVectorToBlock(kSpace, block_data);
  
  // Compute an orthonormal basis from the n_low vectors.
#pragma omp parallel for collapse(2)
  for(int by=0; by<block_scheme[1]; by++) {
    for(int bx=0; bx<block_scheme[0]; bx++) {
      int blk_idx = by * block_scheme[0] + bx;
      
      // Location of the start of the desired block (x runs fastest)
      int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));
      
      for(int i=0; i<n_low; i++) {
	
	// Copy data from the block data into the block ortho array
	for(int k=0; k<block_size; k++) {
	  block_data_ortho[blk_idx][block_size * i + k] = block_data[blk_idx][block_size * i + k];
	}
	
	// Loop up to i
	for(int j=0; j<i; j++) {
	  
	  // Copy data from the jth ortho array into a temp array
	  std::vector<Complex> temp(block_size, 0.0);
	  for(int k=0; k<block_size; k++) temp[k] = block_data_ortho[blk_idx][block_size * j + k];
	  
	  // <Vj|Vi> : inner product
	  Complex ip = 0.0;	  
	  for(int k=0; k<block_size; k++) ip += conj(temp[k]) * block_data[blk_idx][block_size * i + k];
	  
	  // |Vi> = |Vi> - <Vj|Vi>|Vj> : CAXPY project and write to block_data_ortho 
	  for(int k=0; k<block_size; k++) {
	    block_data_ortho[blk_idx][block_size * i + k] = (block_data_ortho[blk_idx][block_size * i + k] - ip * temp[k]);
	  }	  
	}

	// Normalize
	Complex ip = 0.0;
	for(int k=0; k<block_size; k++) {
	  ip += conj(block_data_ortho[blk_idx][block_size * i + k]) * block_data_ortho[blk_idx][block_size * i + k];
	}
	
	double nrm = 1.0/sqrt(ip.real());
	for(int k=0; k<block_size; k++) {
	  block_data_ortho[blk_idx][block_size * i + k] *= nrm;
	}
      }
    }
  }
  
  // Get coefficients: project the blocked n_conv eigenvectors on to the
  // orthonormalised low blocks.
  // Modified Gramm-Schmidt
#pragma omp parallel for collapse(3)
  for(int j=0; j<n_conv; j++) {
    for(int by=0; by<block_scheme[1]; by++) {
      for(int bx=0; bx<block_scheme[0]; bx++) {
	
	int blk_idx = by * block_scheme[0] + bx;        
	for(int i=0; i<n_low; i++) {

	  // Inner product between orthed i block (low) and j block (high)
	  Complex ip = 0.0;
	  for(int k=0; k<block_size; k++) {
	    ip += (conj(block_data_ortho[blk_idx][block_size * i + k]) * 
		   block_data[blk_idx][block_size * j + k]);
	  }
	  
	  block_coeffs[blk_idx][n_low * j + i] = ip;
	  
	  for(int k=0; k<block_size; k++) {
	    block_data[blk_idx][block_size * j + k] = block_data[blk_idx][block_size * j + k] - (ip * block_data_ortho[blk_idx][block_size * i + k]);
	  }
	}
      }
    }
  }
}

void Eig::blockExpand(std::vector<field<Complex>*> &kSpace) {
  
  // Loop over the desired eigenvalues.
#pragma omp parallel for 
  for(int j=0; j<n_conv; j++) {
    
    std::vector<std::vector<Complex>> block_data_temp(n_blocks, std::vector<Complex> (block_size, 0.0));
    
    // Loop over blocks, Gramm-Schmidt the blocks on the low modes
    for(int by=0; by<block_scheme[1]; by++) {
      for(int bx=0; bx<block_scheme[0]; bx++) {
	int blk_idx = by * block_scheme[0] + bx;

	for(int i=0; i<n_low; i++) {	  
	  for(int k=0; k<block_size; k++) {
	    block_data_temp[blk_idx][k] += block_coeffs[blk_idx][j*n_low + i] * block_data_ortho[blk_idx][block_size * i + k];
	  }
	}	
      }
    }

    for(int by=0; by<block_scheme[1]; by++) {
      for(int bx=0; bx<block_scheme[0]; bx++) {
	int blk_idx = by * block_scheme[0] + bx;

	// Location of the start of the desired block (x runs fastest)
	int blk_offset = 2 * ((x_block_size * bx) + Nx * (y_block_size * by));

	// Copy data to the jth eigenvector
	for(int ny=0; ny<y_block_size; ny++) {
	  for(int nx=0; nx<x_block_size; nx++) {
	    for(int mu=0; mu<2; mu++) {
	      
	      int loc_idx = 2*(nx + x_block_size * ny) + mu;      // Local index
	      int glo_idx = blk_offset + 2*((ny * Nx) + nx) + mu; // Global index
	      kSpace[j]->data[glo_idx] = block_data_temp[blk_idx][loc_idx];
	    }
	  }
	}
      }
    }
  }
}
