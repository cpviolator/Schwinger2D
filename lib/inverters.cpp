#include "inverters.h"

inverterCG::inverterCG(param_t param) {
  res = new field<Complex>(param);
  p = new field<Complex>(param);
  Ap = new field<Complex>(param);
  temp = new field<Complex>(param);
  op = MdagM;
}


void inverterCG::OPERATOR(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){
  switch(op) {
  case M: Dpsi(out, in, gauge); break;
  case Mdag: Ddagpsi(out, in, gauge); break;
  case MdagM: DdagDpsi(out, in, gauge); break;
  case MMdag: DDdagpsi(out, in, gauge); break;
  default: cout << "Undefined operator type requested:" << op << endl;
    exit(0);
  }
}

// Get size of shifts from x
int inverterCG::solveMulti(std::vector<field<Complex> *> &x,
			   field<Complex> *b,
			   field<Complex> *gauge, std::vector<double> shifts) {

  double mass_light = gauge->p.m;
  gauge->p.m = gauge->p.m_heavy;
  if(x.size() != shifts.size()) {
    printf("Error: x.size() = %lu, shifts.size() = %lu\n",
	   x.size(), shifts.size());
    exit(0);
  }

  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);
  if(verbose) cout << "bsqrt = " << bsqrt << endl;
  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Warning in inverterCG: inverting on zero source or Nan." << endl;
    exit(0);
  }
  
  int n_shifts = x.size();
  int n_shifts_remaining = n_shifts;
  std::vector<double> zeta(n_shifts, 1.0);
  std::vector<double> zeta_old(n_shifts, 1.0);
  std::vector<double> beta_arr(n_shifts, 0.0);
  std::vector<double> alpha_arr(n_shifts, 1.0);
  std::vector<bool> active(n_shifts, true);
  
  // compute initial residual
  //---------------------------------------  
  res->copy(b);
  p->copy(b);
  rsq = blas::norm2(res->data);
  if(verbose) printf("CG iter 0, rsq = %g\n", rsq);
  
  std::vector<field<Complex> *> p_arr(n_shifts);
  for(int i=0; i<n_shifts; i++) {    
    blas::zero(x[i]->data);
    p_arr[i] = new field<Complex>(x[0]->p);
    p_arr[i]->copy(p);
  }
  
  alpha = 1.0;
  
  // Iterate until convergence
  //---------------------------------------
  bool convergence = false;
  for (iter = 0; iter < gauge->p.max_iter_cg && !convergence; iter++) {
    
    // Compute Ap.
    OPERATOR(Ap, p, gauge);
    denom = (blas::cDotProd(p->data, Ap->data)).real();
    double alpha_old = alpha;
    alpha = -rsq/denom;
    
    // Compute new alpha and zeta
    for (int i=0; i<n_shifts_remaining && active[i]; i++) {

      double c0 = zeta[i] * zeta_old[i] * alpha_old; 
      double c1 = beta * alpha * (zeta_old[i] - zeta[i]);
      double c2 = zeta_old[i] * alpha_old * (1.0 - shifts[i]*alpha);
      
      zeta_old[i] = zeta[i];
      if(c1 + c2 != 0.0) zeta[i] = c0 / (c1 + c2);
      else             zeta[i] = 0.0;
      
      if(zeta[i] != 0.0) alpha_arr[i] = alpha * zeta[i] / zeta_old[i];
      else               alpha_arr[i] = 0.0;
      
      blas::axpy(-alpha_arr[i], p_arr[i]->data, x[i]->data);
    }
    
    blas::axpy(alpha, Ap->data, res->data);    

    rsq_new = blas::norm2(res->data);
    if(verbose) printf("CG iter %d, rsq = %g\n", iter+1, rsq_new);

    // Exit if new residual on the lightest vector is small enough
    bool converged = true;
    for (int i=0; i<n_shifts_remaining && active[i]; i++) {
      if(zeta[i] * sqrt(rsq_new) < gauge->p.eps*bnorm) {
	active[i] = false;
	if(verbose) printf("CG iter %d, shift %d converged: res = %g\n", iter+1, i, zeta[i] * sqrt(rsq_new));
      } else {
	converged = false;
      }
    }
    
    if(converged) {
      rsq = rsq_new;
      convergence = true;
    }
    
    // Update vec using new residual
    beta = rsq_new/rsq;
    rsq = rsq_new;
    
    for (int i=0; i<n_shifts_remaining && active[i]; i++) {
      beta_arr[i] = beta * zeta[i] * alpha_arr[i]/(zeta_old[i] * alpha);
      blas::axpby(zeta[i], res->data, beta_arr[i], p_arr[i]->data);
    }
    
    blas::axpby(1.0, res->data, beta, p->data);
    
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

  if(verbose) {
    //Sanity
    cout << "source norm = " << blas::norm(b->data) << endl;
    for(int i=0; i<n_shifts; i++) {
      cout << "sol norm = " << blas::norm(x[i]->data) << endl;    
      OPERATOR(temp, x[i], gauge);
      blas::caxpy(shifts[i], x[i]->data, temp->data);
      
      blas::axpy(-1.0, temp->data, b->data, res->data);
      double truersq = blas::norm2(res->data);
      printf("CG: Converged shift = %d, rsq = %.16e, truersq = %.16e\n", i, sqrt(rsq), sqrt(truersq)/(bsqrt));
    }
  }
  
  //printf("CG: Converged iter = %d\n", success);
  gauge->p.m = mass_light;
  return success;  
}

int inverterCG::solve(field<Complex> *x, field<Complex> *b,
		      std::vector<field<Complex> *> &kSpace,
		      std::vector<Complex> &evals,
		      field<Complex> *gauge, bool deflate) {
  
  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Warning in inverterCG: inverting on zero source or Nan!" << endl;
    //exit(0);
  }

  blas::zero(temp->data);
  blas::zero(res->data);
  
  // compute initial residual
  //---------------------------------------  
  if (blas::norm2(x->data) > 0.0) {    

    // initial guess supplied: res = b - A*x0
    use_init_guess = true;
    OPERATOR(temp, x, gauge);
    blas::axpy(-1.0, temp->data, b->data, res->data);
    
    // Update bnorm
    bnorm = blas::norm2(res->data);
    bsqrt = sqrt(bnorm);

    // temp contains original guess
    temp->copy(x);
    
    if(verbose) {
      cout << "using initial guess, |x0| = " << blas::norm(temp->data)
	   << ", |b| = " << bsqrt
	   << ", |res| = " << blas::norm(res->data) << endl;
    }
  } else {
    // no initial guess. Initial residual is the source.    
    res->copy(b);
    blas::zero(temp->data);
  }
  
  if (deflate) {
    deflateResidual(temp, res, kSpace, evals);
    OPERATOR(res, temp, gauge);
    blas::axpy(-1.0, res->data, b->data, res->data);
  }
  
  blas::zero(x->data);
  p->copy(res);  
  rsq = blas::norm2(res->data);
  //---------------------------------------
  
  // Iterate until convergence
  //---------------------------------------
  bool convergence = false;
  for (iter = 0; iter < gauge->p.max_iter_cg && !convergence; iter++) {
    
    // Compute Ap.
    OPERATOR(Ap, p, gauge);
 
    denom = (blas::cDotProd(p->data, Ap->data)).real();
    alpha = rsq/denom;
    
    blas::axpy( alpha, p->data,    x->data);
    blas::axpy(-alpha, Ap->data, res->data);
    
    // Exit if new residual is small enough
    rsq_new = blas::norm2(res->data);
    if (verbose) printf("CG iter %d, rsq = %g\n", iter+1, rsq_new);
    //printf("CG iter %d, rsq = %g\n", iter+1, rsq_new);
    if (rsq_new < gauge->p.eps*bnorm) {
      rsq = rsq_new;
      convergence = true;
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
    
    OPERATOR(temp, x, gauge);
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
    cout << "Warning in inverterCG: inverting on zero source or Nan!" << endl;
    exit(0);
  }

  blas::zero(temp->data);
  blas::zero(res->data);
  
  // compute initial residual
  //---------------------------------------  
  if (blas::norm2(x->data) > 0.0) {    

    // initial guess supplied: res = b - A*x0
    use_init_guess = true;
    OPERATOR(temp, x, gauge);
    blas::axpy(-offset, x->data, temp->data);
    blas::axpy(-1.0, temp->data, b->data, res->data);
    
    // Update bnorm
    bnorm = blas::norm2(res->data);
    bsqrt = sqrt(bnorm);

    // temp contains original guess
    temp->copy(x);
    
    if(verbose) {
      cout << "using initial guess, |x0| = " << blas::norm(temp->data)
	   << ", |b| = " << bsqrt
	   << ", |res| = " << blas::norm(res->data) << endl;
    }
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
    OPERATOR(Ap, p, gauge);
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
    
    OPERATOR(temp, x, gauge);
    blas::axpy(-offset, x->data, temp->data);    
    blas::axpy(-1.0, temp->data, b->data, res->data);
    double truersq = blas::norm2(res->data);
    printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", iter+1, rsq, truersq/(bsqrt*bsqrt));
  }
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
