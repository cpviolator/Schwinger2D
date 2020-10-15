#include "inverters.h"

field<Complex> *guess_exp;
field<Complex> *sol1;
field<Complex> *sol2;
field<Complex> *delta_sol;

bool init = false;
int counter = 0;
//bool spline = true;
bool spline = false;

int Ainvpsi(field<Complex> *x, field<Complex> *b, field<Complex> *x0, field<Complex> *gauge) {
  
  int success = 0;

  if(counter > 1 && spline) {
    x0->copy(guess_exp);
  }
  
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

  if (blas::norm2(x0->data) > 0.0) {    
    // initial guess supplied: res = b - A*x0
    
    //Solve the deflated system.
    use_init_guess = true;
    DdagDpsi(temp, x0, gauge);    
    blas::axpy(-1.0, temp->data, res->data);

    // Update bnorm
    //bnorm = blas::norm2(res->data);
    //bsqrt = sqrt(bnorm);
    
    //cout << "using initial guess, |x0| = " << blas::norm(x0->data)
    //<< ", |b| = " << bsqrt
    //<< ", |res| = " << blas::norm(res->data) << endl;
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
    //printf("CG iter %d, rsq = %g\n", k+1, rsq_new);
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
  /*
  DdagDpsi(temp, x, gauge);
  blas::axpy(-1.0, temp->data, b->data, res->data);
  double truersq = blas::norm2(res->data);
  printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", k+1, rsq, truersq/(bsqrt*bsqrt));
  */

  // Sanity
  //cout << "Guess norm = " << blas::norm(x0->data) << endl;
  //cout << "source norm = " << blas::norm(b->data) << endl;
  //cout << "sol norm = " << blas::norm(x->data) << endl;
  
  if(init == false && spline) {
    guess_exp = new field<Complex>(gauge->p);
    sol1 = new field<Complex>(gauge->p);
    sol2 = new field<Complex>(gauge->p);
    delta_sol = new field<Complex>(gauge->p);
    init = true;
  }

  if(counter == 0 && spline) {
    sol1->copy(x);
  }

  if(counter > 0 && spline) {
    sol2->copy(x);
    blas::zero(delta_sol->data);
    // delta_sol = sol2 - sol1 
    blas::axpy( 1.0, sol2->data, delta_sol->data);
    blas::axpy(-1.0, sol1->data, delta_sol->data);
    // sol1 updated
    sol1->copy(sol2);
    // guess updated
    guess_exp->copy(sol2);
    blas::axpy(1.0, delta_sol->data, guess_exp->data);
  }
  
  counter++;
  
  delete res;
  delete p;
  delete Ap;
  delete temp;
  
  return success;  
}

// let dD \equiv (d/dtheta D)
//
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
//
// *****  Should optimize this to operate only on EVEN sites. ****

void forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi, field<Complex> *guess){
  
  if(gauge->p.dynamic == true) {

    blas::zero(fD->data);
    
    //phip = (D^dagD)^-1 * phi
    field<Complex> *phip = new field<Complex>(gauge->p);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1) phi = (D^-1 * g3 * D^-1 g3) phi.
    Ainvpsi(phip, phi, guess, gauge);
    //phip->print();
    
    //g3Dphi = g3D * phip
    field<Complex> *g3Dphi = new field<Complex>(gauge->p); 
    g3Dpsi(g3Dphi, phip, gauge);
    

    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
#pragma omp parallel for
    for(int x=0; x<Nx; x++) {
      int xp1 = (x+1)%Nx;
      //int xm1 = (x-1+Nx)%Nx;
      for(int y=0; y<Ny; y++) {
	int yp1 = (y+1)%Ny;
	//int ym1 = (y-1+Ny)%Ny;	

	double temp = 0.0;
	//mu = 0
	//upper
	// | r  1 | 
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |
	
	temp = real(I*((conj(gauge->read(x,y,0)) *
			 (conj(phip->read(xp1,y,0)) * (r*g3Dphi->read(x,y,0) +   g3Dphi->read(x,y,1)) -
			  conj(phip->read(xp1,y,1)) * (  g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			-
			(gauge->read(x,y,0) *
			 (conj(phip->read(x,y,0)) * (r*g3Dphi->read(xp1,y,0) -   g3Dphi->read(xp1,y,1)) +
			  conj(phip->read(x,y,1)) * (  g3Dphi->read(xp1,y,0) - r*g3Dphi->read(xp1,y,1))))
			)
		     );
	
	fD->write(x,y,0,temp);
	
	//mu = 1
	//upper
	// | r -i | 
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	temp = real(I*((conj(gauge->read(x,y,1)) *
			 (conj(phip->read(x,yp1,0)) * (r*g3Dphi->read(x,y,0) - I*g3Dphi->read(x,y,1)) -
			  conj(phip->read(x,yp1,1)) * (I*g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			-			       
			(gauge->read(x,y,1) *
			 (conj(phip->read(x,y,0)) * (r*g3Dphi->read(x,yp1,0) + I*g3Dphi->read(x,yp1,1)) +
			  conj(phip->read(x,y,1)) * (I*g3Dphi->read(x,yp1,0) - r*g3Dphi->read(x,yp1,1))))
			)
		     );
	
	fD->write(x,y,1,temp);
      }
    }
    
    delete phip;
    delete g3Dphi;    
  }
}
