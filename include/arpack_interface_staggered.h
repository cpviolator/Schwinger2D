/*=====================================================================================

  Thu Aug 02 14:15:21 EDT Dean Howarth
  
  ARPACK interafce for 2D compact U(1) theory.

  ====================================================================================*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>

using namespace std;

void arpackErrorHelpNAUPD();
void arpackErrorHelpNEUPD();

static void mergeAbs(double *sort1, int *idx1, int n1, double *sort2,
		     int *idx2, int n2, bool inverse) {
  int i1=0, i2=0;
  int *ord;
  double *result;
  
  ord    = (int *)    malloc(sizeof(int)   *(n1+n2)); 
  result = (double *) malloc(sizeof(double)*(n1+n2)); 
  
  for(int i=0; i<(n1+n2); i++) {
    if((fabs(sort1[i1]) >= fabs(sort2[i2])) != inverse) { //LOGICAL XOR
      result[i] = sort1[i1];
      ord[i] = idx1[i1];
      i1++;
    } else {
      result[i] = sort2[i2];
      ord[i] = idx2[i2];
      i2++;
    }
    
    if(i1 == n1) {
      for(int j=i+1; j<(n1+n2); j++,i2++) {
	result[j] = sort2[i2];
	ord[j] = idx2[i2];
      }
      i = n1+n2;
    } else if (i2 == n2) {
      for(int j=i+1; j<(n1+n2); j++,i1++) {
	result[j] = sort1[i1];
	ord[j] = idx1[i1];
      }
      i = i1+i2;
    }
  }  
  for(int i=0;i<n1;i++) {
    idx1[i] = ord[i];
    sort1[i] = result[i];
  }
  
  for(int i=0;i<n2;i++) {
    idx2[i] = ord[i+n1];
    sort2[i] = result[i+n1];
  }  
  free (ord);
  free (result);
}
  
static void sortAbs(double *unsorted, int n, bool inverse, int *idx) {
  
  if (n <= 1)
    return;
  
  int n1,n2;
  
  n1 = n>>1;
  n2 = n-n1;
  
  double *unsort1 = unsorted;
  double *unsort2 = (double *)((char*)unsorted + n1*sizeof(double));
  int *idx1 = idx;
  int *idx2 = (int *)((char*)idx + n1*sizeof(int));
  
  sortAbs(unsort1, n1, inverse, idx1);
  sortAbs(unsort2, n2, inverse, idx2);
  
  mergeAbs(unsort1, idx1, n1, unsort2, idx2, n2, inverse);
}



#define ARPACK(s) s ## _

#ifdef __cplusplus
extern "C" {
#endif

  /**
   *  Interface functions to the external ARPACK library. These functions utilize 
   *  ARPACK's implemntation of the Implicitly Restarted Arnoldi Method to compute a 
   *  number of eigenvectors/eigenvalues with user specified features, such as those 
   *  with small real part, small magnitude etc. Parallel (OMP/MPI) versions
   *  are also supported.
   */
  
  
  //Serial, single prec complex eigenvectors
  extern int ARPACK(cnaupd) (int *ido, char *bmat, int *n, char *which, int *nev,
			     float *tol, std::complex<float> *resid, int *ncv,
			     std::complex<float> *v, int *ldv, int *iparam, int *ipntr,
			     std::complex<float> *workd, std::complex<float> *workl,
			     int *lworkl, float *rwork, int *info, int bmat_size,
			     int spectrum_size);
  
  //Serial, double prec complex eigenvectors
  extern int ARPACK(znaupd)(int *ido, char *bmat, int *n, char *which, int *nev,
			    double *tol, std::complex<double> *resid, int *ncv,
			    std::complex<double> *v, int *ldv, int *iparam, int *ipntr,
			    std::complex<double> *workd, std::complex<double> *workl, 
			    int *lworkl, double *rwork, int *info, int bmat_size,
			    int spectrum_size);
  
  //Serial, single prec complex eigenvalues
  extern int ARPACK(cneupd) (int *comp_evecs, char *howmany, int *select,
			     std::complex<float> *evals, std::complex<float> *v,
			     int *ldv, std::complex<float> *sigma,
			     std::complex<float> *workev, char *bmat, int *n,
			     char *which, int *nev, float *tol,
			     std::complex<float> *resid, int *ncv,
			     std::complex<float> *v1, int *ldv1, int *iparam,
			     int *ipntr, std::complex<float> *workd,
			     std::complex<float> *workl, int *lworkl,
			     float *rwork, int *info, int howmany_size, int bmat_size,
			     int spectrum_size);			
  
  //Serial, double prec complex eigenvalues
  extern int ARPACK(zneupd) (int *comp_evecs, char *howmany, int *select,
			     std::complex<double> *evals, std::complex<double> *v,
			     int *ldv, std::complex<double> *sigma,
			     std::complex<double> *workev, char *bmat, int *n,
			     char *which, int *nev, double *tol,
			     std::complex<double> *resid, int *ncv,
			     std::complex<double> *v1, int *ldv1, int *iparam,
			     int *ipntr, std::complex<double> *workd,
			     std::complex<double> *workl, int *lworkl,
			     double *rwork, int *info, int howmany_size, int bmat_size,
			     int spectrum_size);
  
  extern int ARPACK(mcinitdebug)(int*,int*,int*,int*,int*,int*,int*,int*);
    
  //ARPACK initlog and finilog routines for printing the ARPACK log  
  extern int ARPACK(initlog) (int*, char*, int);
  extern int ARPACK(finilog) (int*);
  
#ifdef __cplusplus
}
#endif


void arpackErrorHelpNAUPD() {
  printf("\nError help NAUPD\n\n");
  printf("INFO Integer.  (INPUT/OUTPUT)\n");
  printf("     If INFO .EQ. 0, a randomly initial residual vector is used.\n");
  printf("     If INFO .NE. 0, RESID contains the initial residual vector,\n");
  printf("                        possibly from a previous run.\n");
  printf("     Error flag on output.\n");
  printf("     =  0: Normal exit.\n");
  printf("     =  1: Maximum number of iterations taken.\n");
  printf("        All possible eigenvalues of OP has been found. IPARAM(5)\n");
  printf("        returns the number of wanted converged Ritz values.\n");
  printf("     =  2: No longer an informational error. Deprecated starting\n");
  printf("        with release 2 of ARPACK.\n");
  printf("     =  3: No shifts could be applied during a cycle of the\n");
  printf("        Implicitly restarted Arnoldi iteration. One possibility\n");
  printf("        is to increase the size of NCV relative to NEV.\n");
  printf("        See remark 4 below.\n");
  printf("     = -1: N must be positive.\n");
  printf("     = -2: NEV must be positive.\n");
  printf("     = -3: NCV-NEV >= 1 and less than or equal to N.\n");
  printf("     = -4: The maximum number of Arnoldi update iteration\n");
  printf("        must be greater than zero.\n");
  printf("     = -5: WHICH must be 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
  printf("     = -6: BMAT must be one of 'I' or 'G'.\n");
  printf("     = -7: Length of private work array is not sufficient.\n");
  printf("     = -8: Error return from LAPACK eigenvalue calculation;\n");
  printf("     = -9: Starting vector is zero.\n");
  printf("     = -10: IPARAM(7) must be 1,2,3.\n");
  printf("     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
  printf("     = -12: IPARAM(1) must be equal to 0 or 1.\n");
  printf("     = -9999: Could not build an Arnoldi factorization.\n");
  printf("        User input error highly likely.  Please\n");
  printf("        check actual array dimensions and layout.\n");
  printf("        IPARAM(5) returns the size of the current Arnoldi\n");
  printf("        factorization.\n");
}

void arpackErrorHelpNEUPD() {
  printf("\nError help NEUPD\n\n");
  printf("INFO Integer.  (OUTPUT)\n");
  printf("     Error flag on output.\n");
  printf("     =  0: Normal exit.\n");
  printf("     =  1: The Schur form computed by LAPACK routine csheqr\n");
  printf("        could not be reordered by LAPACK routine ztrsen.\n");
  printf("        Re-enter subroutine zneupd with IPARAM(5)=NCV and\n");
  printf("        increase the size of the array D to have\n");
  printf("        dimension at least dimension NCV and allocate at\n");
  printf("        least NCV\n");
  printf("        columns for Z. NOTE: Not necessary if Z and V share\n");
  printf("        the same space. Please notify the authors if this\n");
  printf("        error occurs.\n");
  printf("     = -1: N must be positive.\n");
  printf("     = -2: NEV must be positive.\n");
  printf("     = -3: NCV-NEV >= 1 and less than or equal to N.\n");
  printf("     = -5: WHICH must be 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
  printf("     = -6: BMAT must be one of 'I' or 'G'.\n");
  printf("     = -7: Length of private work WORKL array is inufficient.\n");
  printf("     = -8: Error return from LAPACK eigenvalue calculation.\n");
  printf("        This should never happened.\n");
  printf("     = -9: Error return from calculation of eigenvectors.\n");
  printf("        Informational error from LAPACK routine ztrevc.\n");
  printf("     = -10: IPARAM(7) must be 1,2,3\n");
  printf("     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
  printf("     = -12: HOWMNY = 'S' not yet implemented\n");
  printf("     = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
  printf("     = -14: ZNAUPD did not find any eigenvalues to sufficient\n");
  printf("        accuracy.\n");
  printf("     = -15: ZNEUPD got a different count of the number of\n");
  printf("        converged Ritz values than ZNAUPD got. This\n");
  printf("        indicates the user probably made an error in\n");
  printf("        passing data from ZNAUPD to ZNEUPD or that the\n");
  printf("        data was modified before entering ZNEUPD\n");
}

/*
void chebyOp(Complex out[LX][LY], Complex in[LX][LY],
	     Complex gauge[LX][LY][2], param_t p ) {
  
  double delta,theta;
  double sigma,sigma1,sigma_old;
  double d1,d2,d3;
  
  double a = 1.0;
  double b = 11.0;
  
  delta = (b-a)/2.0;
  theta = (b+a)/2.0;
  
  sigma1 = -delta/theta;
  
  copyField(out,in);
  
  //if(arpack_param->polyDeg == 0)
  //return;
  
  d1 =  sigma1/delta;
  d2 =  1.0;

  DdagDpsi(out, in, gauge, p);   
  
  axpby(d2, in, d1, out, out);
  
  ////C_1(x) = x
  //if(arpack_param->polyDeg == 1 )
  //return;
  
  // C_0 is the current 'in'  vector.
  // C_1 is the current 'out' vector.
  
  Complex tmp1[LX][LY];
  Complex tmp2[LX][LY];
  
  copyField(tmp1,in);
  copyField(tmp2,out);
  
  //Using Chebyshev polynomial recursion relation,
  //C_{m+1}(x) = 2*x*C_{m} - C_{m-1}
  
  sigma_old = sigma1;
  
  //construct C_{m+1}(x)
  for(int i=2; i<100; i++){
    
    sigma = 1.0/(2.0/sigma1-sigma_old);
    
    d1 = 2.0*sigma/delta;
    d2 = -d1*theta;
    d3 = -sigma*sigma_old;
    
    //mat*C_m
    DdagDpsi(out, in, gauge, p);   
    
    ax(d3,tmp1);
    xpaypbz(tmp1,d2,tmp2,d1,out);
    
    copyField(tmp1,tmp2);
    copyField(tmp2,out);
    sigma_old = sigma;
    
  }
}
*/

int arpack_solve_double(Complex gauge[LX][LY][2], param_t p, Complex guess[LX][LY], int infoGuess, int step, int iter) {
  
  //Construct parameters and memory allocation
  //------------------------------------------
  
  // all FORTRAN communication uses underscored 
  int ido_; 
  int info_;
  int *ipntr_ = (int*)malloc(14*sizeof(int));
  int *iparam_ = (int*)malloc(11*sizeof(int));
  int n_    = LX*LY,
    nev_    = p.nEv,
    nkv_    = p.nKv,
    ldv_    = LX*LY,
    lworkl_ = (3 * nkv_*nkv_ + 5*nkv_) * 2,
    rvec_   = 1;
  int max_iter = p.arpackMaxiter;

  double tol_ = p.arpackTol;

  Complex Zero(0.0,0.0);
  
  double *mod_evals_sorted  = (double*)malloc(nkv_*sizeof(double));
  int *evals_sorted_idx = (int*)malloc(nkv_*sizeof(int));
  
  //Memory checks
  if((mod_evals_sorted == nullptr) ||
     (evals_sorted_idx == nullptr) ) {
    printf("eigenSolver: not enough memory for host eigenvalue sorting");
    exit(0);
  }
  
  //ARPACK workspace
  std::complex<double> sigma_ = 0.0;
  std::complex<double> *resid_ =
    (std::complex<double> *) malloc(ldv_*sizeof(std::complex<double>));
  std::complex<double> *w_workd_ =
    (std::complex<double> *) malloc(3*ldv_*sizeof(std::complex<double>));
  std::complex<double> *w_workl_ =
    (std::complex<double> *) malloc(lworkl_*sizeof(std::complex<double>)); 
  std::complex<double> *w_workev_=
    (std::complex<double> *) malloc(2*nkv_*sizeof(std::complex<double>));    
  double *w_rwork_  = (double *)malloc(nkv_*sizeof(double));    
  int *select_ = (int*)malloc(nkv_*sizeof(int));
  
  std::complex<double> *evecs = (std::complex<double> *) malloc(nkv_*n_*sizeof(std::complex<double>));
  std::complex<double> *evals = (std::complex<double> *) malloc(nkv_   *sizeof(std::complex<double>));

  for(int n=0; n<nkv_; n++) {
    evals[n] = Zero;
    for(int x=0; x<LX; x++) {
      for(int y=0; y<LY; y++) {
	evecs[n*LX*LY + x*LY + y] = Zero;
	if(n==0) resid_[y*LX + x] = guess[x][y];
      }
    }
  }

  //for(int k=0; k<100; k++) cout << resid_[k] << " ";
  //cout << endl;
  
  //Alias pointers
  std::complex<double> *evecs_ = nullptr;
  evecs_ = (std::complex<double>*) (double*)(evecs);    
  std::complex<double> *evals_ = nullptr;
  evals_ = (std::complex<double>*) (double*)(evals);
  
  //Memory checks
  if((iparam_ == nullptr) ||
     (ipntr_ == nullptr) || 
     (resid_ == nullptr) ||  
     (w_workd_ == nullptr) || 
     (w_workl_ == nullptr) ||
     (w_workev_ == nullptr) ||
     (w_rwork_ == nullptr) || 
     (select_ == nullptr) ) {
    printf("eigenSolver: not enough memory for ARPACK workspace.\n");
    exit(0);
  }    

  //Assign values to ARPACK params 
  ido_        = 0;
  info_       = infoGuess;
  iparam_[0]  = 1;
  iparam_[1]  = 1;
  iparam_[2]  = max_iter;
  iparam_[3]  = 1;
  iparam_[6]  = 1;
  //iparam_[7]  = 1;
  
  //ARPACK problem type to be solved
  char howmny='P';
  char bmat = 'I';
  char *spectrum;
  spectrum = strdup("SR"); //Initialsed just to stop the compiler warning...
  if(p.polyACC) spectrum = strdup("LR");
  int iter_cnt= 0;

  //Start ARPACK routines
  //---------------------------------------------------------------------------------
 
  Complex *psi1;
  Complex *psi2;

  Complex psi1_cpy[LX][LY];
  Complex psi2_cpy[LX][LY];
  
  for(int x=0; x<LX; x++) {
    for(int y=0; y<LY; y++) {
      psi1_cpy[x][y] = Zero;
      psi2_cpy[x][y] = Zero;
    }
  }
  
  psi1 = w_workd_;
  psi2 = w_workd_ + n_;

  double t1;
  double time = 0.0;;
  do {
    
    t1 = -((double)clock());
    
    //Interface to arpack routines
    //----------------------------
    
    ARPACK(znaupd)(&ido_, &bmat, &n_, spectrum, &nev_, &tol_, resid_, &nkv_,
		   evecs_, &n_, iparam_, ipntr_, w_workd_, w_workl_, &lworkl_,
		   w_rwork_, &info_, 1, 2);

    if (info_ != 0) {
      printf("\nError in znaupd info = %d. Exiting...\n",info_);
      arpackErrorHelpNAUPD();
      exit(0);
    }
    
    if (ido_ == 99 || info_ == 1)
      break;
    
    if (ido_ == -1 || ido_ == 1) {

      //Copy from Arpack workspace
      for(int x=0; x<LX; x++) {
	for(int y=0; y<LY; y++) {
	  psi1_cpy[x][y] = *(psi1 + y*LX + x);
	}
      }
      //Apply matrix-vector operation
      //if(p.polyACC == 1) chebyOp(psi2_cpy, psi1_cpy, gauge, p);
      DdagDpsi(psi2_cpy, psi1_cpy, gauge, p);
      
      //Copy to Arpack workspace
      for(int x=0; x<LX; x++) {
	for(int y=0; y<LY; y++) {
	  *(psi2 + y*LX + x) = psi2_cpy[x][y];
	}
      }
    }
    
    t1 += clock();
    time += t1;
    if(iter_cnt % 100 == 0) printf("Arpack Iteration: %d (%e secs)\n", iter_cnt, time/CLOCKS_PER_SEC);
    iter_cnt++;
    
  } while (99 != ido_ && iter_cnt < max_iter);
  
  //Subspace calulated sucessfully. Compute nEv eigenvectors and values   
  printf("Finished in %e secs: iter=%04d  info=%d  ido=%d\n", time/CLOCKS_PER_SEC, iter_cnt, info_, ido_);      
  //printf("Computing eigenvectors\n");
  if(infoGuess == 1) {
    for(int x=0; x<LX; x++) {
      for(int y=0; y<LY; y++) {
	guess[x][y] = resid_[y*LX + x];
      }
    }
  }
  
  //Interface to arpack routines
  //----------------------------
  ARPACK(zneupd)(&rvec_, &howmny, select_, evals_, evecs_, &n_, &sigma_,
		 w_workev_, &bmat, &n_, spectrum, &nev_, &tol_,
		 resid_, &nkv_, evecs_, &n_, iparam_, ipntr_, w_workd_,
		 w_workl_, &lworkl_, w_rwork_, &info_, 1, 1, 2);
  if (info_ == -15) {
    printf("\nError in zneupd info = %d. You likely need to\n"
	   "increase the maximum ARPACK iterations. Exiting...\n", info_);
    arpackErrorHelpNEUPD();
    exit(0);
  } else if (info_ != 0) {
    printf("\nError in zneupd info = %d. Exiting...\n", info_);
    arpackErrorHelpNEUPD();
  }
  
  int nconv = iparam_[4];
  for(int j=0; j<nconv; j++){    
    evals_sorted_idx[j] = j;
    mod_evals_sorted[j] = std::abs(evals_[j]);
  }
  
  //Sort the Ritz in absolute ascending order
  t1 =  -((double)clock());
  bool inverse = (p.polyACC == 1 ? false : true);
  sortAbs(mod_evals_sorted, nconv, inverse, evals_sorted_idx);
  t1 +=  clock();
  
  printf("Sorting time: %f sec\n",t1/CLOCKS_PER_SEC);
  printf("Sorted eigenvalues based on their absolute values:\n");
  
  // Print additional convergence information.
  if( (info_) == 1){
    printf("Maximum number of iterations reached.\n");
  }
  else{
    if(info_ == 3){
      printf("Error: No shifts could be applied during implicit\n");
      printf("Error: Arnoldi update, try increasing NkV.\n");
    }
  }
  
  
  //Print Evalues
  t1 = -(double)clock();
  Complex psi3[LX][LY];
  
  for(int i=0; i<nev_ ;i++){    
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) {
	psi3[x][y] = evecs[evals_sorted_idx[i]*n_ + y*LX + x];
      }
    
    //apply matrix-vector operation here:
    DdagDpsi(psi2_cpy, psi3, gauge, p);
    
    // lambda = v^dag * M*v    
    evals_[i] = dotField(psi3, psi2_cpy);
    
    Complex unit(1.0,0.0);
    Complex m_lambda(-real(evals_[i]),
		     -imag(evals_[i]));
    
    // d_v = ||M*v - lambda*v||
    caxpby(unit, psi2_cpy, m_lambda, psi3, psi1_cpy);

    double L2norm = norm2(psi1_cpy);    
    printf("Eval[%04d]: %+.16e %+.16e Residual: %+.3e\n",
	   i, real(evals_[i]), imag(evals_[i]), sqrt(L2norm));    
  }    
  
  t1 += clock();
  printf("Eigenvalues of Dirac operator computed in: %f sec\n", t1/CLOCKS_PER_SEC);
  
  // cleanup 
  free(ipntr_);
  free(iparam_);
  free(mod_evals_sorted);
  free(evals_sorted_idx);
  free(resid_);
  free(w_workd_);
  free(w_workl_);
  free(w_workev_);
  free(w_rwork_);
  free(select_);
  free(spectrum);
  
  return iter_cnt;  
}  
