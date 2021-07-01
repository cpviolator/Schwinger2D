#include "inverters.h"
#include "iram.h"

#define OPERATOR DdagDpsi

inverterCG::inverterCG(param_t param) {
  res = new field<Complex>(param);
  p = new field<Complex>(param);
  Ap = new field<Complex>(param);
  temp = new field<Complex>(param);    
}

int inverterCG::solve(field<Complex> *x, field<Complex> *b,
		      std::vector<field<Complex> *> &kSpace, std::vector<Complex> &evals,
		      field<Complex> *gauge, bool deflate) {
  
  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Warning in inverterCG: inverting on zero source or (Nan!" << endl;
    //exit(0);
  }

  blas::zero(temp->data);
  blas::zero(res->data);
  
  // compute initial residual
  //---------------------------------------  
  if (blas::norm2(x->data) > 0.0) {    

    // initial guess supplied: res = b - A*x0
    use_init_guess = true;
    DdagDpsi(temp, x, gauge);    
    blas::axpy(-1.0, temp->data, b->data, res->data);
    
    // Update bnorm
    bnorm = blas::norm2(res->data);
    bsqrt = sqrt(bnorm);

    // temp contains original guess
    temp->copy(x);
    
    //if(verbose) {
    //cout << "using initial guess, |x0| = " << blas::norm(temp->data)
    //<< ", |b| = " << bsqrt
    //<< ", |res| = " << blas::norm(res->data) << endl;
    //}
  } else {
    // no initial guess. Initial residual is the source.    
    res->copy(b);
    blas::zero(temp->data);
  }
  
  if (deflate) {
    deflateResidual(temp, res, kSpace, evals);
    DdagDpsi(res, temp, gauge);
    blas::axpy(-1.0, res->data, b->data, res->data);
  }
  
  blas::zero(x->data);
  p->copy(res);  
  rsq = blas::norm2(res->data);
  //---------------------------------------
  
  // Iterate until convergence
  //---------------------------------------
  for (iter = 0; iter < gauge->p.max_iter_cg; iter++) {
    
    // Compute Ap.
    DdagDpsi(Ap, p, gauge);
    
    denom = (blas::cDotProd(p->data, Ap->data)).real();
    alpha = rsq/denom;
    
    blas::axpy( alpha, p->data,    x->data);
    blas::axpy(-alpha, Ap->data, res->data);
    
    // Exit if new residual is small enough
    rsq_new = blas::norm2(res->data);
    if (verbose) printf("CG iter %d, rsq = %g\n", iter+1, rsq_new);
    if (rsq_new < gauge->p.eps*bnorm) {
      rsq = rsq_new;
      break;
    }
    
    // Update vec using new residual
    beta = rsq_new/rsq;
    rsq = rsq_new;
    
    blas::axpy(beta, p->data, res->data, p->data);
    
  } // End loop over iter
  //---------------------------------------

  if(iter == gauge->p.max_iter_cg) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", iter+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = iter+1; 
  }

  // x contains the solution to the deflated system b - A*x0.
  // We must add back the exact part if using an initial guess too
  blas::axpy(1.0, temp->data, x->data);  
  // x now contains the solution to the RHS b.
  
  if(verbose) {
    
    //Sanity
    cout << "source norm = " << blas::norm(b->data) << endl;
    cout << "sol norm = " << blas::norm(x->data) << endl;
    
    DdagDpsi(temp, x, gauge);
    blas::axpy(-1.0, temp->data, b->data, res->data);
    double truersq = blas::norm2(res->data);
    printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", iter+1, rsq, truersq/(bsqrt*bsqrt));
  }
  //printf("CG: Converged iter = %d\n", success);  
  return success;  
}

int inverterCG::solve(field<Complex> *x, field<Complex> *b, field<Complex> *gauge) {
  std::vector<field<Complex> *> kSpace;
  std::vector<Complex> evals;
  return solve(x, b, kSpace, evals, gauge, false);
}

int inverterCG::solve(field<Complex> *x, field<Complex> *b, field<Complex> *gauge, double offset) {

  std::vector<field<Complex> *> kSpace;
  std::vector<Complex> evals;

  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Warning in inverterCG: inverting on zero source or (Nan!" << endl;
    //exit(0);
  }

  blas::zero(temp->data);
  blas::zero(res->data);
  
  // compute initial residual
  //---------------------------------------  
  if (blas::norm2(x->data) > 0.0) {    

    // initial guess supplied: res = b - A*x0
    use_init_guess = true;
    DdagDpsi(temp, x, gauge);
    blas::axpy(-offset, x->data, temp->data);
    
    blas::axpy(-1.0, temp->data, b->data, res->data);
    
    // Update bnorm
    bnorm = blas::norm2(res->data);
    bsqrt = sqrt(bnorm);

    // temp contains original guess
    temp->copy(x);
    
    //if(verbose) {
    //cout << "using initial guess, |x0| = " << blas::norm(temp->data)
    //<< ", |b| = " << bsqrt
    //<< ", |res| = " << blas::norm(res->data) << endl;
    //}
  } else {
    // no initial guess. Initial residual is the source.    
    res->copy(b);
    blas::zero(temp->data);
  }
  
  blas::zero(x->data);
  p->copy(res);  
  rsq = blas::norm2(res->data);
  //---------------------------------------
  
  // Iterate until convergence
  //---------------------------------------
  for (iter = 0; iter < gauge->p.max_iter_cg; iter++) {
    
    // Compute Ap.
    DdagDpsi(Ap, p, gauge);
    blas::axpy(-offset, p->data, Ap->data);
    
    denom = (blas::cDotProd(p->data, Ap->data)).real();
    alpha = rsq/denom;
    
    blas::axpy( alpha, p->data,    x->data);
    blas::axpy(-alpha, Ap->data, res->data);
    
    // Exit if new residual is small enough
    rsq_new = blas::norm2(res->data);
    if (verbose) printf("CG iter %d, rsq = %g\n", iter+1, rsq_new);
    if (rsq_new < gauge->p.eps*bnorm) {
      rsq = rsq_new;
      break;
    }
    
    // Update vec using new residual
    beta = rsq_new/rsq;
    rsq = rsq_new;
    
    blas::axpy(beta, p->data, res->data, p->data);
    
  } // End loop over iter
  //---------------------------------------

  if(iter == gauge->p.max_iter_cg) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", iter+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = iter+1; 
  }

  // x contains the solution to the deflated system b - A*x0.
  // We must add back the exact part if using an initial guess too
  blas::axpy(1.0, temp->data, x->data);  
  // x now contains the solution to the RHS b.
  
  if(verbose) {
    
    //Sanity
    cout << "source norm = " << blas::norm(b->data) << endl;
    cout << "sol norm = " << blas::norm(x->data) << endl;
    
    DdagDpsi(temp, x, gauge);
    blas::axpy(-offset, x->data, temp->data);
    
    blas::axpy(-1.0, temp->data, b->data, res->data);
    double truersq = blas::norm2(res->data);
    printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", iter+1, rsq, truersq/(bsqrt*bsqrt));
  }
  //printf("CG: Converged iter = %d\n", success);  
  return success;  
}

void inverterCG::deflateResidual(field<Complex> *phi_defl, field<Complex> *phi,
				 std::vector<field<Complex> *> &kSpace,
				 std::vector<Complex> &evals) {
  
  Complex scalar;
  //Deflate each converged eigenpair from the phi
  // phi_defl = (v * lambda^-1 * v^dag) * phi
  for(int i=0; i<kSpace.size(); i++) {
    //Compute scalar part: s = (lambda)^-1 * (v^dag * phi)
    scalar = blas::cDotProd(kSpace[i]->data, phi->data);
    scalar /= evals[i].real();
    //Accumulate in phi_defl: phi_defl = phi_defl + v * s
    blas::caxpy(scalar, kSpace[i]->data, phi_defl->data);
  }
}

int Ainvpsi(field<Complex> *x, field<Complex> *b, field<Complex> *x0, field<Complex> *gauge) {
  
  int success = 0;
  bool verbose = false;
  
  field<Complex> *res = new field<Complex>(gauge->p);
  field<Complex> *p = new field<Complex>(gauge->p);
  field<Complex> *Ap = new field<Complex>(gauge->p);
  field<Complex> *temp = new field<Complex>(gauge->p);
  
  double alpha = 0.0, beta = 0.0, denom = 0.0, rsq = 0.0, rsq_new = 0.0, bsqrt = 0.0, bnorm = 0.0;
  bool use_init_guess = false;  
  
  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Error in Wilson Ainvpsi: inverting on zero source... or nan!" << endl;
    exit(0);
  }

  // compute initial residual 
  res->copy(b);
  
  // initial guess supplied: res = b - A*x0
  if (blas::norm2(x0->data) > 0.0) {    
    
    use_init_guess = true;
    DdagDpsi(temp, x0, gauge);    
    blas::axpy(-1.0, temp->data, res->data);

    // Update bnorm
    bnorm = blas::norm2(res->data);
    bsqrt = sqrt(bnorm);

    if(verbose) {
      cout << "using initial guess, |x0| = " << blas::norm(x0->data)
	   << ", |b| = " << bsqrt
	   << ", |res| = " << blas::norm(res->data) << endl;
    }
  }
  
  p->copy(res);
  
  rsq = blas::norm2(res->data);
  
  // Iterate until convergence
  int k;
  for (k=0; k<gauge->p.max_iter_cg; k++) {
    
    // Compute Ap.
    DdagDpsi(Ap, p, gauge);
    
    denom = (blas::cDotProd(p->data, Ap->data)).real();
    alpha = rsq/denom;
    
    blas::axpy( alpha, p->data,    x->data);
    blas::axpy(-alpha, Ap->data, res->data);
    
    // Exit if new residual is small enough
    rsq_new = blas::norm2(res->data);
    if (verbose) printf("CG iter %d, rsq = %g\n", k+1, rsq_new);
    if (rsq_new < gauge->p.eps*bnorm) {
      rsq = rsq_new;
      break;
    }
    
    // Update vec using new residual
    beta = rsq_new/rsq;
    rsq = rsq_new;
    
    blas::axpy(beta, p->data, res->data, p->data);
    
  } // End loop over k
  
  if(k == gauge->p.max_iter_cg) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", k+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = k+1; 
  }
  
  if(use_init_guess) {
    // x contains the solution to the deflated system b - A*x0.
    // We must add back the exact part
    blas::axpy(1.0, x0->data, x->data);
    // x now contains the solution to the RHS b.
  }
  
  //printf("CG: Converged iter = %d\n", k+1);
  
  if(verbose) {
    
    // Sanity
    cout << "Guess norm = " << blas::norm(x0->data) << endl;
    cout << "source norm = " << blas::norm(b->data) << endl;
    cout << "sol norm = " << blas::norm(x->data) << endl;
    
    DdagDpsi(temp, x, gauge);
    blas::axpy(-1.0, temp->data, b->data, res->data);
    double truersq = blas::norm2(res->data);
    printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", k+1, rsq, truersq/(bsqrt*bsqrt));
  }

    
  delete res;
  delete p;
  delete Ap;
  delete temp;
  
  return success;  
}
