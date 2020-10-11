#include "utils.h"
#include "inverters.h"

int Ainvpsi(field<Complex> *x, field<Complex> *b, field<Complex> *x0, field<Complex> *gauge) {
  
  int success = 0;

  field<Complex> *res = new field<Complex>(gauge->p);
  field<Complex> *p = new field<Complex>(gauge->p);
  field<Complex> *Ap = new field<Complex>(gauge->p);
  field<Complex> *temp = new field<Complex>(gauge->p);
  
  double alpha = 0.0, beta = 0.0, denom = 0.0, rsq = 0.0, rsq_new = 0.0, bsqrt = 0.0, bnorm = 0.0;
  bool deflating = false;
  
  // Find norm of rhs.
  bnorm = blas::norm2(b->data);
  bsqrt = sqrt(bnorm);
  
  if(bsqrt == 0 || bsqrt != bsqrt) {
    cout << "Error in Wilson Ainvpsi: inverting on zero source... or nan!" << endl;
    exit(0);
  }
  res->copy(b);
  
  // res = b - A*x0
  if (blas::norm2(x0->data) != 0.0) {
    
    //Solve the deflated system.
    deflating = true;
    DdagDpsi(temp, x0, gauge);    
    blas::caxpy(-1.0, temp->data, res->data);
    
    cout << "using initial guess, |x0| = " << blas::norm(x0->data)
	 << ", |b| = " << bsqrt
	 << ", |res| = " << blas::norm(res->data) << endl;
  }
  
  p->copy(res);
  rsq = blas::norm(res->data);
  
  // Iterate until convergence
  int k;
  for (k=0; k<gauge->p.max_iter_cg; k++) {
    
    // Compute Ap.
    DdagDpsi(Ap, p, gauge);
    
    denom = blas::dotProd(p->data, Ap->data).real();
    alpha = rsq/denom;
    
    blas::caxpy( alpha, p->data,    x->data);
    blas::caxpy(-alpha, Ap->data, res->data);
    
    // Exit if new residual is small enough
    rsq_new = blas::norm2(res->data);
    //printf("CG iter %d, rsq = %g\n", k+1, rsq_new);
    if (rsq_new < gauge->p.eps*bnorm) {
      rsq = rsq_new;
      break;
    }
    
    // Update vec using new residual
    beta = rsq_new/rsq;
    rsq = rsq_new;
    
    blas::caxpy(beta, p->data, res->data, p->data);
    
  } // End loop over k
  
  if(k == gauge->p.max_iter_cg) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", k+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = 1; 
  }
  
  if(deflating) {
    // x contains the solution to the deflated system b - A*x0.
    // We must add back the exact part
    blas::caxpy(1.0, x0->data, x->data);
    // x now contains the solution to the RHS b.
  }

  //delete res;
  //delete p;
  //delete Ap;
  //delete temp;
  
  //DdagDpsi(temp, x, gauge);
  //blas::caxpy(-1.0, temp->data, b, res->data);
  //double truersq = real(dotField(res, res));
  //printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", k+1, rsq, truersq/(bsqrt*bsqrt));
  return success;
  
}

// let dD \equiv (d/dtheta D)
//
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
//
// *****  Should optimize this to operate only on EVEN sites. ****

void forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi, field<Complex> *guess){
  
  if(gauge->p.dynamic == true) {

    //blas::zero(fD->data);
    
    //phip = (D^dagD)^-1 * phi
    field<Complex> *phip = new field<Complex>(gauge->p);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1) phi = (D^-1 * g3 * D^-1 g3) phi.
    field<Complex> *guess = new field<Complex>(gauge->p);
    Ainvpsi(phip, phi, guess, gauge);
    
    //g3Dphi = g3D * phip
    field<Complex> *g3Dphi = new field<Complex>(gauge->p); 
    g3Dpsi(g3Dphi, phip, gauge);
    
    double temp = 0.0;
    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
    for(int x=0; x<Nx; x++) {
      for(int y=0; y<Ny; y++) {

	//mu = 0
	//upper
	// | r  1 | 
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |
	
	temp = real(I*((conj(gauge->read(x,y,0)) *
			 (conj(phip->read(x+1,y,0)) * (r*g3Dphi->read(x,y,0) +   g3Dphi->read(x,y,1)) -
			  conj(phip->read(x+1,y,1)) * (  g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			-
			(gauge->read(x,y,0) *
			 (conj(phip->read(x,y,0)) * (r*g3Dphi->read(x+1,y,0) -   g3Dphi->read(x+1,y,1)) +
			  conj(phip->read(x,y,1)) * (  g3Dphi->read(x+1,y,0) - r*g3Dphi->read(x+1,y,1))))
			)
		     );
	
	fD->write(x,y,1,temp);
	
	//mu = 1
	//upper
	// | r -i | 
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	temp = real(I*((conj(gauge->read(x,y,1)) *
			 (conj(phip->read(x,y+1,0)) * (r*g3Dphi->read(x,y,0) - I*g3Dphi->read(x,y,1)) -
			  conj(phip->read(x,y+1,1)) * (I*g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			-			       
			(gauge->read(x,y,1) *
			 (conj(phip->read(x,y,0)) * (r*g3Dphi->read(x,y+1,0) + I*g3Dphi->read(x,y+1,1)) +
			  conj(phip->read(x,y,1)) * (I*g3Dphi->read(x,y+1,0) - r*g3Dphi->read(x,y+1,1))))
			)
		     );
	
	fD->write(x,y,1,temp);
      }
    }
    
    //delete phip;
    //delete guess;
    //delete g3Dphi;    
  }
}
