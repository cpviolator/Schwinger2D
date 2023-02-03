#include "iram.h"

IRAM::IRAM(EigParam eig_param_in) {
  
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
  iram_verbose = eig_param.iram_verbose;
  op = eig_param.op;
  block_scheme[0] = eig_param.block_scheme[0];
  block_scheme[1] = eig_param.block_scheme[1];
  n_low = eig_param.n_low;

  x_block_size = Nx/block_scheme[0];
  y_block_size = Ny/block_scheme[1];
  n_blocks = block_scheme[0]*block_scheme[1];
  block_size = 2 * x_block_size * y_block_size;

  if(Nx%block_scheme[0] != 0) {
    cout << "IRAM: Error: x block_scheme = " << block_scheme[0] << " does not divide Nx = " << Nx << endl;
    exit(0);
  }
  if(Ny%block_scheme[1] != 0) {
    cout << "IRAM: Error: y block_scheme = " << block_scheme[1] << " does not divide Ny = " << Ny << endl;
    exit(0);
  }

  // n_blocks * n_conv * (block_size + l_low) complex elems
  for(int i=0; i<n_blocks; i++) {
    block_data_ortho.push_back(std::vector<Complex> (n_low * block_size, 0.0));
    block_coeffs.push_back(std::vector<Complex> (n_low * n_conv, 0.0));
  }
}

void IRAM::OPERATOR(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){
  switch(op) {
  case M: Dpsi(out, in, gauge); break;
  case Mdag: Ddagpsi(out, in, gauge); break;
  case MdagM: DdagDpsi(out, in, gauge); break;
  case MMdag: DdagDpsi(out, in, gauge); break;
  default: cout << "IRAM: Undefined operator type requested: " << op << endl;
    exit(0);
  }
}

void IRAM::deflate(field<Complex> *guess, field<Complex> *phi) {

  if(kSpace_defl.size() == 0) {
    cout << "IRAM Error: No deflation space evecs exist" << endl;
    exit(0);
  }
  if(evals_defl.size() == 0) {
    cout << "IRAM: Error: No deflation space evals exist" << endl;
    exit(0);
  }
  
  blas::zero(guess->data);
  Complex scalar;
  
  // Deflate each converged eigenpair from the guess  
  // guess_defl = (v * lambda^-1 * v^dag) * guess
  if(verbosity) cout << "IRAM: Deflating " << n_deflate << " eigenpairs." << endl;
  for(int i=0; i<n_deflate; i++) {
    
    if(use_compressed_space) { 
      // Compute scalar part: s = (lambda)^-1 * (v^dag * phi)
      scalar = blas::cDotProd(kSpace_mg[i]->data, phi->data);
      scalar /= real(evals_mg[i]);
      // Accumulate in guess defl_guess: defl_guess = defl_guess + v * s
      blas::caxpy(scalar, kSpace_mg[i]->data, guess->data);
    }
    else {
      // Compute scalar part: s = (lambda)^-1 * (v^dag * phi)
      scalar = blas::cDotProd(kSpace_defl[i]->data, phi->data);
      scalar /= real(evals_defl[i]);
      // Accumulate in guess defl_guess: defl_guess = defl_guess + v * s
      blas::caxpy(scalar, kSpace_defl[i]->data, guess->data);
    }
  }
}

void IRAM::zsortc(Spectrum which, int n, std::vector<Complex> &x, std::vector<Complex> &y) {
  
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
  default: cout << "IRAM: Undefined sort" << endl;
  }

  // Repopulate x and y arrays with sorted elements
  for(int i=0; i<n; i++) {
    x[i] = array[i].first;
    y[i] = array[i].second;
  }
}

// Overloaded version of zsortc to deal with real y array.
void IRAM::zsortc(Spectrum which, int n, std::vector<Complex> &x, std::vector<double> &y) {

  std::vector<Complex> y_tmp(n,0.0);
  for(int i=0; i<n; i++) y_tmp[i].real(y[i]);
  zsortc(which, n, x, y_tmp);
  for(int i=0; i<n; i++) y[i] = y_tmp[i].real();
}

// Overloaded version of zsortc to deal with real x array.
void IRAM::zsortc(Spectrum which, int n, std::vector<double> &x, std::vector<Complex> &y) {

  std::vector<Complex> x_tmp(n,0.0);
  for(int i=0; i<n; i++) x_tmp[i].real(x[i]);
  zsortc(which, n, x_tmp, y);
  for(int i=0; i<n; i++) x[i] = x_tmp[i].real();
}

// Overloaded version of zsortc to deal with real x and y array.
void IRAM::zsortc(Spectrum which, int n, std::vector<double> &x, std::vector<double> &y) {

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



void IRAM::arnoldiStep(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
		       Eigen::MatrixXcd &upperHessEigen,
		       field<Complex> *r, double &beta, int j) {

  //%---------------------------------------------------%
  //| STEP 1: Check if the B norm of j-th residual      |
  //| vector is zero. Equivalent to determine whether   |
  //| an exact j-step Arnoldi factorization is present. |
  //%---------------------------------------------------%
  beta = blas::norm(r->data);
  //for(int i=0; i<10; i++) cout << (r.data)[i] << endl;
  //exit(0);
  //beta = blas::norm2(r.data.data(), r.data.size());
  //cout << beta << endl;
  
  
  //%--------------------------------%
  //| STEP 2:  v_{j} = r_{j-1}/rnorm |
  //%--------------------------------%  
  blas::ax(1.0/beta, r->data);
  blas::copy(kSpace[j]->data, r->data);
  
  //%----------------------------%
  //| STEP 3:  r_{j} = OP*v_{j}; |
  //%----------------------------%  
  OPERATOR(r, kSpace[j], gauge);
  
  //%-------------------------------------%
  //| The following is needed for STEP 5. |
  //| Compute the B-norm of OP*v_{j}.     |
  //%-------------------------------------%

  double wnorm = blas::norm(r->data);

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
    upperHessEigen(i,j) = blas::cDotProd(kSpace[i]->data, r->data);
  }
  
  //%--------------------------------------%
  //| Orthogonalize r_{j} against V_{j}.   |
  //| RESID contains OP*v_{j}. See STEP 3. | 
  //%--------------------------------------%
  //r = r - H_{j,i} * v_j 
  for (int i = 0; i < j+1; i++) {
    blas::caxpy(-1.0*upperHessEigen(i,j), kSpace[i]->data, r->data);
  }
  
  if(j > 0) upperHessEigen(j,j-1) = beta;
  
  beta = blas::norm(r->data);
  
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
      alpha[i] = blas::cDotProd(kSpace[i]->data, r->data);
      upperHessEigen(i,j) += alpha[i];
    }
    
    for(int i=0; i < j+1; i++) {
      blas::caxpy(-1.0*alpha[i], kSpace[i]->data, r->data);
    }
    
    beta = blas::norm(r->data);
    orth_iter++;
  }
  
  if(orth_iter == orth_iter_max) {
    //%---------------------------------------%
    //| RESID is numerically in the span of V |
    //%---------------------------------------%
    cout << "IRAM: Unable to orthonormalise r" << endl;
    exit(0);
  }
}

void IRAM::reorder(std::vector<field<Complex> *> kSpace, std::vector<Complex> evals, std::vector<double> residua, int n_kr, Spectrum spectrum) {

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
  default: printf("Undefined spectrum type %d given", spectrum);
    exit(0);
  }
  
  // Repopulate arrays with sorted elements
  for(int i=0; i<n; i++) {
    std::swap(evals[i], std::get<0>(array[i]));
    std::swap(residua[i], std::get<1>(array[i]));
    std::swap(kSpace[i], std::get<2>(array[i]));
  }
}

void IRAM::eigensolveFromUpperHess(MatrixXcd &upperHessEigen, MatrixXcd &Qmat,
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

void IRAM::qriteration(MatrixXcd &Rmat, MatrixXcd &Qmat, const int n_kr, const double tol)
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

int IRAM::qrFromUpperHess(MatrixXcd &upperHess, MatrixXcd &Qmat, std::vector<Complex> &evals, std::vector<double> &residua, const double beta, const int n_kr, const double tol)
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
  
  if(iram_verbose) printf("IRAM: eigensystem iterations = %d\n", iter);

  return iter;  
}

void IRAM::rotateVecsComplex(std::vector<field<Complex> *> vecs, Eigen::MatrixXcd mat, int num_locked, int iter_keep, int n_kr) {

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

void IRAM::computeEvals(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace, std::vector<double> &residua, std::vector<Complex> &evals, int n_ev) {
  
  //temp vector
  field<Complex> *temp = new field<Complex>(gauge->p);
  for (int i = 0; i < n_ev; i++) {
    
    // r = A * v_i
    OPERATOR(temp, kSpace[i], gauge);

    // lambda_i = v_i^dag A v_i / (v_i^dag * v_i)
    evals[i] = blas::cDotProd(kSpace[i]->data, temp->data);

    // Measure ||lambda_i*v_i - A*v_i||
    Complex n_unit(-1.0, 0.0);
    blas::caxpby(evals[i], kSpace[i]->data, n_unit, temp->data);
    residua[i] = blas::norm(temp->data);
  }
  delete temp;
}

void IRAM::iram(const field<Complex> *gauge, std::vector<field<Complex> *> kSpace,
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

  // START init
  //---------------------------------------------------------  
  gettimeofday(&start, NULL);  

  if (!(n_kr > n_ev + 6) && n_kr != (gauge->p.Nx * gauge->p.Ny * 2)) {
    printf("IRAM WARNING: n_kr=%d must be greater than n_ev+6=%d\n", n_kr, n_ev + 6);
    exit(0);
  }

  if(iram_verbose) {
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
  bool converged = false;
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
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
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
  r->copy(kSpace[0]);
  
  // START IRAM
  // Implicitly restarted Arnoldi method for asymmetric eigenvalue problems
  if(iram_verbose) printf("IRAM: Start IRAM solution\n");
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
	if(iram_verbose) printf("IRAM: residuum[%d] = %e, cond = %e\n", i, residua[idx], tol * abs(evals[idx]));
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
      cout << "No shifts can be applied" << endl;
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

      blas::caxpby(upperHessEigen(num_keep, num_keep-1), kSpace[num_keep]->data, Qmat(n_kr-1, num_keep-1), r->data);
      gettimeofday(&end, NULL);  
      t_compute += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

      if(blas::norm(r->data) < epsilon) {
	printf("IRAM: Congratulations! You have reached an invariant subspace at iter %d, beta = %e\n", restart_iter, blas::norm(r->data));
	exit(0);
      }
    }
    restart_iter++;    
  }

  kSpace.resize(n_conv);
  evals.resize(n_conv);
  
  gettimeofday(&total_end, NULL);  
  t_total = ((total_end.tv_sec - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
  
  // Post computation report  
  if (!converged) {    
    printf("IRAM: Failed to compute the requested %d vectors with with a %d search space and a %d Krylov space in %d restart_steps and %d OPs.\n", n_conv, n_ev, n_kr, restart_iter, iter);
    exit(0);
  } else {
    if(verbosity) {
      printf("IRAM: Computed the requested %d vectors with a %d search space and a %d Krylov space in %d restart_steps and %d OPs in %e secs.\n", n_conv, n_ev, n_kr, restart_iter, iter, (t_compute + t_sort + t_EV + t_QR));    
      for (int i = 0; i < n_conv; i++) {
	printf("IRAM: EigValue[%04d]: ||(%+.8e, %+.8e)|| = %+.8e residual %.8e\n", i, evals[i].real(), evals[i].imag(), abs(evals[i]), residua[i]);
      }
    }
  }

  if(iram_verbose) {
    cout << "IRAM: Timings:" << endl;
    cout << "IRAM: init = " << t_init << endl;
    cout << "IRAM: compute = " << t_compute << endl;
    cout << "IRAM: sort = " << t_sort << endl;
    cout << "IRAM: EV = " << t_EV << endl;
    cout << "IRAM: QR = " << t_QR << endl;
    cout << "IRAM: missing = " << (t_total) << " - " << (t_compute + t_init + t_sort + t_EV + t_QR + t_eigen) << " = " << (t_total - (t_compute + t_init + t_sort + t_EV + t_QR + t_eigen)) << " ("<<(100*((t_total - (t_compute + t_init + t_sort + t_EV + t_QR + t_eigen))))/t_total<<"%)" << endl;  
  }
}

void IRAM::inspectEvolvedSpectrum(const field<Complex> *gauge, int iter) {
  
  std::vector<field<Complex>*> kSpace;
  std::vector<field<Complex>*> kSpace_aux;
  std::vector<Complex> evals;

  // Allocate space to hold the previous kSpace
  if(kSpace_pre.size() == 0) prepareKrylovSpace(kSpace_pre, evals_pre, gauge->p);  

  // Allocate space for this inspection
  prepareKrylovSpace(kSpace, evals, gauge->p);

  // Compute eigenspectrum of this gauge.
  if(inspection_counter == 0) {
    // A deflation space was just constructed, copy it to the
    // kSpace container and save.
    cout << "IRAM: First inspection call, copy deflation space" << endl;
    for(int i=0; i<n_conv; i++) {
      blas::copy(kSpace[i]->data, kSpace_defl[i]->data);
      evals[i] = evals_defl[i];
    }
    
  } else {
    // If inspection_counter > 0, the gauge field has changed and
    // we need to recompute the space for inspecton
    cout << "IRAM: Computing new spectrum" << endl;
    iram(gauge, kSpace, evals);
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
    for(int i=0; i<n_conv; i++) {
      for(int j=0; j<n_conv; j++) inner_products[i + n_conv * j] = blas::cDotProd(kSpace_pre[j]->data, kSpace[i]->data);
    }
    
    // Measure kSpace overlaps
    std::vector<double> overlaps(n_conv, 0.0);
    for(int i=0; i<n_conv; i++) 
      for(int j=0; j<n_conv; j++) overlaps[i] += (inner_products[i + n_conv * j] * conj(inner_products[i + n_conv * j])).real();
    
    // Dump inner product data
    name = "data/eig/inner_products" + to_string(iter);
    constructName(name, gauge->p);
    name += ".dat";
    snprintf(fname, 100, "%s", name.c_str());	  
    fp = fopen(fname, "a");
    for(int i=0; i<n_conv; i++) {
      for(int j=0; j<n_conv; j++) {
	fprintf(fp, "%d   %d %d %.16e %.16e\n",
		inspection_counter, i, j,
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
  
  // Compute an new MG deflation space 
  if(inspection_counter > 0) {
    prepareKrylovSpace(kSpace_aux, evals, gauge->p);
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
  
  // Copy the new space into the old space for overlap comparision.
  for(int i=0; i<n_conv; i++) {
    blas::copy(kSpace_pre[i]->data, kSpace[i]->data);
    evals_pre[i] = evals[i];
  }
  
  // Clean up
  for(int i=0; i<kSpace.size(); i++) delete kSpace[i];
  for(int i=0; i<kSpace_aux.size(); i++) delete kSpace_aux[i];
  
  inspection_counter++;
}

void IRAM::computeDeflationSpace(const field<Complex> *gauge){
  prepareKrylovSpace(kSpace_defl, evals_defl, gauge->p);
  iram(gauge, kSpace_defl, evals_defl);

  prepareKrylovSpace(kSpace_mg, evals_mg, gauge->p);
  computeMGDeflationSpace(kSpace_mg, kSpace_defl, gauge);
}

void IRAM::computeMGDeflationSpace(std::vector<field<Complex>*> &kSpace_out, const std::vector<field<Complex>*> &kSpace_in,
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

  for (int i = 0; i < n_conv; i++) {
    printf("IRAM: Post Compression EigValue[%04d]: ||(%+.8e, %+.8e)|| = %+.8e residual %.8e\n", i, evals_mg[i].real(), evals_mg[i].imag(), abs(evals_mg[i]), resid[i]);
  }
  
  // Check compression ratio
  int pre = 2 * n_conv * Nx * Ny;
  int post = n_blocks * n_low * (block_size + n_conv);
  cout << "IRAM: Algorithmic compression report: " << endl;
  cout << "IRAM: Complex(double) elems pre = " << pre << " Complex(double) elems post = " << post << endl;
  cout << "IRAM: Ratio: " << (100.0 * post)/pre << "% of original data." << endl;
  cout << "IRAM: " << n_low << " low eigenvectors used " << endl;
  cout << "IRAM: " << (n_conv - n_low) << " high eigenvectors reconstructed " << endl;
  cout << "IRAM: Compress/decompress time = " << ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 << endl;
}

void IRAM::prepareKrylovSpace(std::vector<field<Complex>*> &kSpace,
			      std::vector<Complex> &evals,
			      const Param p) {  
  
  kSpace.reserve(eig_param.n_conv);
  for(int i=0; i<eig_param.n_conv; i++) kSpace.push_back(new field<Complex>(p));
  evals.resize(eig_param.n_conv);
}

// Read the block data from the (iEig)th vector in kSpace 
void IRAM::readVectorToBlock(const std::vector<field<Complex>*> &kSpace, std::vector<std::vector<Complex>> &block_data) {
  
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

void IRAM::blockCompress(const std::vector<field<Complex>*> &kSpace) {
  
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

void IRAM::blockExpand(std::vector<field<Complex>*> &kSpace) {
  
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
