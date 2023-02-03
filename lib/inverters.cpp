#include "inverters.h"

inverterCG::inverterCG(Param param) {
  res = new field<Complex>(param);
  p = new field<Complex>(param);
  Ap = new field<Complex>(param);
  temp = new field<Complex>(param);
  op = MdagM;

  verbosity = param.verbosity;
  cg_verbosity = param.cg_verbosity;  
  if(param.current_hmc_iter > param.therm) {    
    deflate = param.eig_param.n_deflate > 0 ? true : false;
    inspect_spectrum = param.inspect_spectrum;
  }
  else deflate = false;
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

  // Switch to heavy mass
  double mass_light = gauge->p.m;
  gauge->p.m = gauge->p.m_heavy;

  // Sanity check
  if(x.size() != shifts.size()) {
    printf("CG Multi: Error: x.size() = %lu, shifts.size() = %lu\n",
	   x.size(), shifts.size());
    exit(0);
  }

  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);
  if(cg_verbosity) cout << "CG Multi: bsqrt = " << bsqrt << endl;
  
  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "CG Multi: Warning in inverterCG Multi: inverting on zero source or Nan." << endl;
    exit(0);
  }
  
  int n_shifts = x.size();
  int n_shifts_remaining = n_shifts;
  std::vector<double> zeta(n_shifts, 1.0);
  std::vector<double> zeta_old(n_shifts, 1.0);
  std::vector<double> beta_arr(n_shifts, 0.0);
  std::vector<double> alpha_arr(n_shifts, 1.0);
  std::vector<bool> active(n_shifts, true);
  
  // Compute initial residual
  //---------------------------------------  
  res->copy(b);
  p->copy(b);
  rsq = blas::norm2(res->data);
  if(cg_verbosity) printf("CG Multi: iter 0, res = %g\n", sqrt(rsq));
  
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
    if(cg_verbosity) printf("CG Multi: iter %d, res = %g\n", iter+1, sqrt(rsq_new));

    // Exit if new residual on the lightest vector is small enough
    bool converged = true;
    for (int i=0; i<n_shifts_remaining && active[i]; i++) {
      if(zeta[i] * sqrt(rsq_new) < pow(gauge->p.tol_cg, 2)*bnorm) {
	active[i] = false;
	if(cg_verbosity) printf("CG Multi: iter %d, shift %d converged: res = %g\n", iter+1, i, zeta[i] * sqrt(rsq_new));
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
    printf("CG Multi: Failed to converge iter = %d, res = %.16e\n", iter, sqrt(rsq)); 
    success = 0; 
  } else {
    // Convergence 
    success = iter; 
  }

  if(cg_verbosity) {
    //Sanity
    cout << "CG Multi: source norm = " << blas::norm(b->data) << endl;
    for(int i=0; i<n_shifts; i++) {
      cout << "CG Multi: sol " << i << " norm = " << blas::norm(x[i]->data) << endl;
    }
    
    for(int i=0; i<n_shifts; i++) {
      OPERATOR(temp, x[i], gauge);
      blas::caxpy(shifts[i], x[i]->data, temp->data);
      
      blas::axpy(-1.0, temp->data, b->data, res->data);
      double truersq = blas::norm2(res->data);
      printf("CG Multi: Converged shift = %d, iter = %d, res = %.16e, trueres = %.16e\n", i, success, sqrt(rsq), sqrt(truersq)/(bsqrt));
    }
  }

  // Restore light mass
  gauge->p.m = mass_light;
  return success;  
}

int inverterCG::solve(field<Complex> *x, const field<Complex> *b, const field<Complex> *gauge) {
  
  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) cout << "CG: Warning in inverterCG: inverting on zero source or Nan!" << endl;
  
  blas::zero(temp->data);
  blas::zero(res->data);

  if(deflate) {
    if (!eig) {
      if(verbosity) cout << "CG: Computing deflation space for CG" << endl;
      eig = new IRAM(gauge->p.eig_param);
      eig->computeDeflationSpace(gauge);
    }
    if(inspect_spectrum) {
      if(verbosity) cout << "CG: Computing spectrum for inspection" << endl;
      eig->inspectEvolvedSpectrum(gauge, gauge->p.current_hmc_iter);
    }
  }
  
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
    
    if(cg_verbosity) {
      cout << "CG: Using initial guess, |x0| = " << blas::norm(temp->data)
	   << ", |b| = " << bsqrt
	   << ", |res| = " << blas::norm(res->data) << endl;
    }
  } else {
    // no initial guess. Initial residual is the source.
    if(cg_verbosity) {
      cout << "CG: No initial guess, |b| = |res| = " << bsqrt;
    }
    res->copy(b);
    blas::zero(temp->data);
  }

  if (deflate && eig->deflationSpaceExists()) {
    if(verbosity) cout << "CG: Applying deflation space preconditioner" << endl;
    eig->deflate(temp, res);
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
    if (cg_verbosity) printf("CG: iter %d, res = %g\n", iter, sqrt(rsq_new));
    if (rsq_new < pow(gauge->p.tol_cg, 2)*bnorm) {
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
    printf("CG: Failed to converge iter = %d, res = %.16e\n", iter, sqrt(rsq)); 
    success = 0; 
  } else {
    // Convergence 
    success = iter; 
  }

  // x contains the solution to the deflated system b - A*x0.
  // We must add back the exact part if using an initial guess too
  blas::axpy(1.0, temp->data, x->data);  
  // x now contains the solution to the RHS b.
  
  if(cg_verbosity) {    
    //Sanity
    cout << "CG: source norm = " << blas::norm(b->data) << endl;
    cout << "CG: sol norm = " << blas::norm(x->data) << endl;
  }    
  OPERATOR(temp, x, gauge);
  blas::axpy(-1.0, temp->data, b->data, res->data);
  double truersq = blas::norm2(res->data);
  if(verbosity) printf("CG: Converged iter = %d, res = %.16e, truersq = %.16e\n", success, sqrt(rsq), sqrt(truersq)/bsqrt);
  
  return success;  
}

int inverterCG::solve(field<Complex> *x, field<Complex> *b, field<Complex> *gauge, double offset) {

  std::vector<field<Complex> *> kSpace;
  std::vector<Complex> evals;

  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);

  // Sanity  
  if(bsqrt == 0 || bsqrt != bsqrt) cout << "CG: Warning in inverterCG: inverting on zero source or Nan!" << endl;

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
    
    if(cg_verbosity) {
      cout << "CG: using initial guess, |x0| = " << blas::norm(temp->data)
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
    if (cg_verbosity) printf("CG: iter %d, rsq = %g\n", iter, rsq_new);
    if (rsq_new < pow(gauge->p.tol_cg, 2)*bnorm) {
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
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", iter, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = iter; 
  }

  // x contains the solution to the deflated system b - A*x0.
  // We must add back the exact part if using an initial guess too
  blas::axpy(1.0, temp->data, x->data);  
  // x now contains the solution to the RHS b.
  
  if(cg_verbosity) {    
    //Sanity
    cout << "CG: source norm = " << blas::norm(b->data) << endl;
    cout << "CG: sol norm = " << blas::norm(x->data) << endl;
  }
  
  OPERATOR(temp, x, gauge);
  blas::axpy(-offset, x->data, temp->data);    
  blas::axpy(-1.0, temp->data, b->data, res->data);
  double truersq = blas::norm2(res->data);
  if(verbosity) printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", iter, rsq, truersq/(bsqrt*bsqrt));

  return success;  
}
