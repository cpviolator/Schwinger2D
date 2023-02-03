#pragma once

#include "schwinger2d_internal.h"
#include "utils.h"
#include "dirac_op.h"
#include "iram.h"
#include "blas.h"

//===============================================================
// CG solutions to Apsi = b 
// see http://en.wikipedia.org/wiki/Conjugate_gradient_method
//===============================================================
//       x: The solution
//       b: The RHS vector
//      x0: An initial guess
//   gauge: The gauge field defining the operator
//   param: The parameter container

// Wilson g5Dg5D matrix inverter
//---------------------------------------------------------------
class inverterCG {

private:
  int success;
  bool verbosity;
  bool cg_verbosity;
  
  field<Complex> *res;
  field<Complex> *p;
  field<Complex> *Ap;
  field<Complex> *temp;
  
  double alpha, beta, denom, rsq, rsq_new, bsqrt, bnorm;
  bool use_init_guess;
  bool deflate;
  bool inspect_spectrum;
  int iter;

  Operator op;

  IRAM *eig;
  
public:

  // Class instance constructor  
  inverterCG(Param param);

  // Operator to solve
  void OPERATOR(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);

  // CG solver
  int solve(field<Complex> *x, const field<Complex> *b, const field<Complex> *gauge);
  
  // CG solver With constant offset 
  int solve(field<Complex> *x, field<Complex> *b, field<Complex> *gauge, double offset);

  // Multi RHS for RHMC
  int solveMulti(std::vector<field<Complex> *> &x, field<Complex> *b,
		 field<Complex> *gauge, std::vector<double> shifts);

  void switchOffDeflation() {
    if(verbosity) cout << "CG: Switching off deflation" << endl;
    deflate = false;
    inspect_spectrum = false;
  };
  void switchOnDeflation() {
    if(verbosity) cout << "CG: Switching on deflation" << endl;
    deflate = true;
    inspect_spectrum = true;
  };
  
  ~inverterCG() {
    delete res;
    delete p;
    delete Ap;
    delete temp;
    if(eig) delete eig;
  }
  
};



/*
//Staggered
int Ainvpsi(Complex psi[LX][LY], Complex b[LX][LY], Complex psi0[LX][LY], const Complex gauge[LX][LY][2], param_t p) {

  int success = 0;
  
  Complex res[LX][LY] , pvec[LX][LY], Apvec[LX][LY];
  double alpha, beta, denom ;
  double rsq = 0, rsqNew = 0, bsqrt = 0.0;
  
  //Intialize  
  zeroField(res);
  zeroField(Apvec);  
  zeroField(pvec);
  
  // Find norm of rhs.
  bsqrt = real(dotField(b,b));
  bsqrt = sqrt(bsqrt);
  
  copyField(res, b); // res = b  - A psi0, for now start with phi0 = 0
  copyField(pvec, res);

  rsq = real(dotField(res,res));
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge, p);

  // iterate till convergence
  int k;
  for (k=0; k<p.maxIterCG; k++) {
    
    denom = real(dotField(pvec,Apvec));
    alpha = rsq/denom;

    axpy( alpha, pvec, psi);
    axpy(-alpha, Apvec, res);
    
    // Exit if new residual is small enough
    rsqNew = real(dotField(res,res));
    
    if (sqrt(rsqNew) < p.eps*bsqrt) {
      //printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // Update vec using new residual
    beta = rsqNew / rsq;
    rsq = rsqNew;
    
    axpy(beta, pvec, res, pvec);
    
    // Compute the new Ap.
    DdagDpsi(Apvec, pvec, gauge,p);  
  }
  //End loop over k

  if(k == p.maxIterCG) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    success = 0; // Failed convergence 
  }
  else {
    success = 1; // Convergence 
    k++;
  }

  DdagDpsi(Apvec, psi, gauge,p);
  axpy(-1.0, Apvec, b, res);
  
  //double truersq =  real(dotField(res,res));
  //printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
}


void forceD(double fD[LX][LY][2], Complex gauge[LX][LY][2], Complex phi[LX][LY],
	    Complex guess[LY][LX], param_t p) {

  if(p.dynamic == true) {

    zeroLat(fD);
    
    Complex phip[LX][LY];
    Complex Dphip[LX][LY];
    zeroField(phip);
    zeroField(Dphip);
    
    Ainvpsi(phip, phi, phip, gauge, p); // note phip = 0 for ODD
    Dpsi(Dphip, phip, gauge, p);        // restrict to Dslash, m = 0
    
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++){
	if( (x+y)%2 == 1) phip[x][y]  = Complex(0.0,0.0);
	if( (x+y)%2 == 0) Dphip[x][y] = Complex(0.0,0.0);
      }
    
    double eta1;
    int yp1, xp1;
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) {
	
	eta1 =(1.0 - 2.0*(x%2));
	xp1 = (x+1)%LX;
	yp1 = (y+1)%LY;
	
	if( (x+y+1)%2 == 0){ 
	  fD[x][y][0] += 2.0*imag(conj(Dphip[x][y]) * gauge[x][y][0] * phip[xp1][y]);
	}
	else {
	  fD[x][y][0] += 2.0*imag(conj(Dphip[xp1][y]) * conj(gauge[x][y][0]) * phip[x][y]);
	};
	
	if( (x+y+1)%2 == 0){    
	  fD[x][y][1] += 2.0*eta1*imag(conj(Dphip[x][y]) * gauge[x][y][1] * phip[x][yp1]);
	}
	else {
	  fD[x][y][1] += 2.0*eta1*imag(conj(Dphip[x][yp1]) * conj(gauge[x][y][1]) * phip[x][y]);
	}
      }	    
  }
}
*/

