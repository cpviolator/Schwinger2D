#include "measurements.h"

void measWilsonLoops(field<Complex> *gauge, double plaq, int iter){

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  
  std::vector<std::vector<Complex>> wLoops(Nx/2 ,std::vector<Complex> (Ny/2 ,0.0));  
  std::vector<double> sigma(Nx/2, 0.0);
  
  Complex w;
  int p1, p2, dx, dy, x, y;
  double inv_Lsq = 1.0/(Nx*Ny);  
  int loop_max = gauge->p.loop_max;

  //Smear the gauge field
  field<Complex> *smeared = new field<Complex>(gauge->p);  
  smearLink(smeared, gauge);
  
  //Loop over all X side sizes of rectangle 
  for(int Xrect=1; Xrect<loop_max; Xrect++) {
      
    //Loop over all Y side sizes of rectangle
    for(int Yrect=1; Yrect<loop_max; Yrect++) {
      
      //Loop over all x,y starting points
      for(x=0; x<Nx; x++)
	for(y=0; y<Ny; y++){
	    
	  w = Complex(1.0,0.0);
	    
	  //Move in +x up to p1.
	  for(dx=0; dx<Xrect; dx++)     w *= smeared->read(x+dx,y,0);
	    
	  //Move in +y up to p2 (p1 constant)
	  p1 = (x + Xrect)%Nx;
	  for(dy=0; dy<Yrect; dy++)     w *= smeared->read(p1,y+dy,1);
	  
	  //Move in -x from p1 to (p2 constant)
	  p2 = (y + Yrect)%Ny;
	  for(dx=Xrect-1; dx>=0; dx--)  w *= conj(smeared->read(x+dx,p2,0));
	  
	  //Move in -y from p2 to y
	  for(dy=Yrect-1; dy>=0; dy--)  w *= conj(smeared->read(x,y+dy,1));
	  wLoops[Xrect][Yrect] += w*inv_Lsq;
	}
    }
  }
  
  //Compute string tension
  for(int size=1; size<loop_max; size++) {
    sigma[size]  = -log(abs((real(wLoops[size][size])/real(wLoops[size-1][size]))* 
			    (real(wLoops[size-1][size-1])/real(wLoops[size][size-1]))));
    
    sigma[size] += -log(abs((real(wLoops[size][size])/real(wLoops[size][size-1]))* 
			    (real(wLoops[size-1][size-1])/real(wLoops[size-1][size]))));
    
    sigma[size] *= 0.5;
    
  }

  string name;
  char fname[256];
  FILE *fp;
  
  name = "data/creutz/creutz";
  constructName(name, gauge->p);
  name += ".dat";
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d %.16e ", iter+1, -log(abs(plaq)) );
  for(int size=2; size<loop_max; size++)
    fprintf(fp, "%.16e ", sigma[size]);
  fprintf(fp, "\n");
  fclose(fp);
  
  for(int sizex=2; sizex<loop_max; sizex++)
    for(int sizey=sizex-1; (sizey < loop_max && sizey <= sizex+1); sizey++) {
      name = "data/rect/rectWL";
      name += "_" + to_string(sizex) + "_" + to_string(sizey);
      constructName(name, gauge->p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %.16e %.16e\n", iter+1, real(wLoops[sizex][sizey]), imag(wLoops[sizex][sizey]));	    
      fclose(fp);
    }
  
  return;
}

/*
//Polyakov loops. x is the spatial dim, y is the temporal dim.
void measPolyakovLoops(Complex gauge[LX][LY][2], int iter, param_t p){

  Complex pLoops[LX/2];
  for(int x=0; x<LX/2; x++) pLoops[x] = 0.0;
  
  Complex w1, w2;
  //Eack polyakov loop correlation is defined by its delta x value.
  //We start at x0, separate to x0 + (x0+L/2-1), and loop over all
  //x0=1 -> x0 = L/2-1.

  //Starting x
  for(int x=0; x<p.loopMax; x++) {

    //Loop over time
    w1 = Complex(1.0,0.0);
    for(int dy=0; dy<LY; dy++) w1 *= gauge[x][dy][1];
    
    //x separation
    for(int dx=0; dx<LX/2; dx++) {
      
      w2 = Complex(1.0,0.0);
      for(int dy=0; dy<LY; dy++) w2 *= gauge[x+dx][dy][1];
      
      pLoops[dx] += conj(w1)*w2/(1.0*LX/2);
      
    }
  }

  string name;
  char fname[256];
  FILE *fp;
  
  name= "data/polyakov/polyakov";
  constructName(name, p);
  name += ".dat";
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  for(int size=1; size<p.loopMax; size++)
    fprintf(fp, "%.16e %.16e ",
	    real(pLoops[size]),
	    imag(pLoops[size]) );
  fprintf(fp, "\n");
  fclose(fp);
  
  name = "data/polyakov/polyakov_ratios";
  constructName(name, p);
  name += ".dat";  
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  for(int size=1 ; size < p.loopMax-1; size++)
    fprintf(fp, "%.16e ",
	    real(pLoops[size+1])/real(pLoops[size]));
  fprintf(fp, "\n");
  fclose(fp);
  
  return;
}

//Pion correlation function
//                              |----------------|
//                              |        |-------|---------|
//  < pi(x) | pi(0) > = < ReTr[dn*(x) g3 up(x) | dn(0) g3 up*(0)] >     
//                    = < ReTr(g3 Gd[0,x] g3 Gu[x,0]) >  
//                    = < ReTr(G*[x,0] G[x,0]) >
//
// using g3 G[x,0] g3 = G*[x,0] and Gu[x,0] \eq Gd[x,0]
//
// if H = Hdag, Tr(H * Hdag) = Sum_{n,m} (H_{n,m}) * (H_{n,m})^*,
// i.e., the sum of the modulus squared of each element

void measPionCorrelation(Complex gauge[LX][LY][2], int top, int iter, param_t p) {
  
  //Up type fermion prop
  Complex propUp[LX][LY][2];
  //Down type fermion prop
  Complex propDn[LX][LY][2];
  //fermion prop CG guess
  Complex propGuess[LX][LY][2];
  //Deflation eigenvectors
  Complex defl_evecs[NEV][LX][LY][2];
  //Deflation eigenvalues
  Complex defl_evals[NEV];
      
  double pion_corr[LY];
            
  Complex source[LX][LY][2];
  Complex Dsource[LX][LY][2];

  char fname[256];
  string name;
  FILE *fp;
  
  //Deflate if requested
#ifdef USE_ARPACK
  arpack_solve(gauge, defl_evecs, defl_evals, 0, 0, p);
#endif
  
  //Up type source
  zeroField(source);
  zeroField(Dsource);
  zeroField(propUp);
  zeroField(propGuess);
  source[0][0][0] = cUnit;
  
  // up -> (g3Dg3) * up *****
  // (g3Dg3D)^-1 * (g3Dg3) up = D^-1 * up *****

  g3psi(source);
  g3Dpsi(Dsource, source, gauge, p);  

  if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
  Ainvpsi(propUp, Dsource, propGuess, gauge, p);

  //Down type source
  zeroField(source);
  zeroField(Dsource);
  zeroField(propDn);
  zeroField(propGuess);
  source[0][0][1] = cUnit;	    
      
  // dn -> (g3Dg3) * dn *****
  // (g3Dg3D)^-1 * (g3Dg3) dn = D^-1 * dn ***** 

  g3psi(source);
  g3Dpsi(Dsource, source, gauge, p);

  if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
  Ainvpsi(propDn, Dsource, propGuess, gauge, p);

  //Get estimate of vacuum trace
  Complex q[2] = {0.0,0.0};
  q[0] = propUp[0][0][0];
  q[1] = propDn[0][0][1];

  double vacuum_trace = (q[0] - q[1]).real();
  
  name = "data/vacuum/estimate_vacuum_Q" + std::to_string(abs(top));
  constructName(name, p);
  name += ".dat";
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  fprintf(fp, "%.16e\n", vacuum_trace);
  fclose(fp);
  
  //Let y be the 'time' dimension
  double corr = 0.0, tmp = 0.0;
  for(int y=0; y<LY; y++) {
    //initialise
    pion_corr[y] = 0.0;
    //Loop over space and spin, fold propagator
    corr = 0.0;
    for(int x=0; x<LX; x++) {
      tmp = abs((conj(propDn[x][y][0]) * propDn[x][y][0] +
      		 conj(propDn[x][y][1]) * propDn[x][y][1] +
      		 conj(propUp[x][y][0]) * propUp[x][y][0] +
      		 conj(propUp[x][y][1]) * propUp[x][y][1]));

      corr += tmp;
    }
    
    //Compute folded propagator
    if ( y < ((LY/2)+1) ) pion_corr[y] += corr;
    else {
      pion_corr[LY-y] += corr;
      pion_corr[LY-y] /= 2.0;
    }
  }
  
  //topological sector pion correlation
  name = "data/pion/pion_Q" + std::to_string(abs(top));
  constructName(name, p);
  name += ".dat";  
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  for(int t=0; t<LY/2+1; t++)
    fprintf(fp, "%.16e ", pion_corr[t]);
  fprintf(fp, "\n");
  fclose(fp);

  //Full pion correlation
  name = "data/pion/pion";
  constructName(name, p);
  name += ".dat";  
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  for(int t=0; t<LY/2+1; t++)
    fprintf(fp, "%.16e ", pion_corr[t]);
  fprintf(fp, "\n");
  fclose(fp);

}

void measVacuumTrace(Complex gauge[LX][LY][2], int top, int iter, param_t p) {
  
  //Up type fermion prop
  Complex propUp[LX][LY][2];
  //Down type fermion prop
  Complex propDn[LX][LY][2];
  //fermion prop CG guess
  Complex propGuess[LX][LY][2];
  //Deflation eigenvectors
  Complex defl_evecs[NEV][LX][LY][2];
  //Deflation eigenvalues
  Complex defl_evals[NEV];
  
  double vacuum_trace[2] = {0.0, 0.0};
  
  Complex source[LX][LY][2];
  Complex Dsource[LX][LY][2];
  
  //Disconnected
  //Loop over time slices
  for(int y=0; y<LY; y++) {
    
    //Sum over spatial sites
    for(int x=0; x<LX; x++) {
      
      //Up type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propUp);
      source[x][y][0] = cUnit;
      
      g3psi(source);
      g3Dpsi(Dsource, source, gauge, p);
      if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
      Ainvpsi(propUp, Dsource, propGuess, gauge, p);
      
      //Down type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propDn);
      source[x][y][1] = cUnit;
      
      g3psi(source);
      g3Dpsi(Dsource, source, gauge, p);
      if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
      Ainvpsi(propDn, Dsource, propGuess, gauge, p);
      
      vacuum_trace[0] += (conj(propDn[x][y][0]) * propDn[x][y][0] +
			  conj(propDn[x][y][1]) * propDn[x][y][1] +
			  conj(propUp[x][y][0]) * propUp[x][y][0] +
			  conj(propUp[x][y][1]) * propUp[x][y][1]).real();
      
      vacuum_trace[1] += (conj(propDn[x][y][0]) * propDn[x][y][0] +
			  conj(propDn[x][y][1]) * propDn[x][y][1] +
			  conj(propUp[x][y][0]) * propUp[x][y][0] +
			  conj(propUp[x][y][1]) * propUp[x][y][1]).imag();
      
    }    
  }

  string name = "data/vacuum/vacuum_Q" + std::to_string(abs(top));
  constructName(name, p);
  name += ".dat";
  
  char fname[256];
  sprintf(fname, "%s", name.c_str());
  FILE *fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  fprintf(fp, "%.16e %.16e\n", vacuum_trace[0], vacuum_trace[1]);
  fclose(fp);
  
}

double measTopCharge(Complex gauge[LX][LY][2], param_t p){
  
  Complex w;
  double top = 0.0;  
  Complex smeared[LX][LY][2];
  smearLink(smeared, gauge, p);
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      w = (smeared[x][y][0] * smeared[ (x+1)%LX ][y][1] *
	   conj(smeared[x][ (y+1)%LY ][0])*conj(smeared[x][y][1]));
      top += arg(w);  // -pi < arg(w) < pi  Geometric value is an integer.
      //print local def here for topology dynamics
      //printf("arg(w) = [ arg(link1) + arg(link2) + c_arg(link3) + c_arg(link4)]
    }
  return top/TWO_PI;
}
*/

double measGaugeAction(field<Complex> *gauge) {

  double beta = gauge->p.beta;
  double action = 0.0;
  Complex plaq = 0.0;

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;  
  for(int x=0; x<Nx;x++)
    for(int y=0; y<Ny; y++){      
      plaq = (gauge->read(x,y,0) * gauge->read(x+1,y,1) *
	      conj(gauge->read(x,y+1,0)) * conj(gauge->read(x,y,1)));
      
      action += beta*real(1.0 - plaq);      
    }  
  return action;
}

double measMomAction(field<double> *mom) {
  
  double action = 0.0;
  double temp = 0.0;
  int Nx = mom->p.Nx;
  int Ny = mom->p.Ny;
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++){
      for(int mu=0; mu<2; mu++){
	temp = mom->read(x,y,mu);
	action += 0.5 * temp * temp;
      }
    }
  
  return action;
}

//Staggered fermion
/*
double measFermAction(Complex gauge[LX][LY][2], Complex phi[LX][LY],
		      param_t p, bool postStep) {
  
  double Hferm = 0.0;
  Complex phitmp[LX][LY];
  
  // cout << "Before Fermion force H = " << H << endl;
  Complex scalar = Complex(0.0,0.0);
  zeroField(phitmp);
  Ainvpsi(phitmp, phi, phitmp, gauge, p);
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      if((x+y)%2 == 0)
	scalar += conj(phi[x][y])*phitmp[x][y];
    }
  
  Hferm += real(scalar);
  //cout << "After Fermion Force H  = " << H << endl;
  
  return Hferm;
}


//Staggered Action
double measAction(double mom[LX][LY][2], Complex gauge[LX][LY][2],
		  Complex phi[LX][LY], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) H += measFermAction(gauge, phi, p, postStep);
  
  return H;
}
*/

//Wilson fermion
double measFermAction(field<Complex> *gauge, field<Complex> *phi, bool postStep) {
  
  field<Complex> *phi_tmp = new field<Complex>(gauge->p);  
  
  //cout << "Before Fermion force H = " << H << endl;
  Complex action = 0.0;
  
  if(postStep) Ainvpsi(phi_tmp, phi, phi_tmp, gauge);
  else phi_tmp->copy(phi);

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++){
      for(int mu=0; mu<2; mu++){
	action += conj(phi->read(x,y,mu))*phi_tmp->read(x,y,mu);
      }
    }    
  
  return action.real();
}

//Wilson Action
double measAction(field<double> *mom, field<Complex> *gauge, field<Complex> *phi, bool postStep) {
  
  double action = 0.0;
  action += measMomAction(mom);
  action += measGaugeAction(gauge);
  if (gauge->p.dynamic) action += measFermAction(gauge, phi, postStep);
  
  return action;
}


Complex measPlaq(field<Complex> *gauge) {
  Complex plaq = 0.0;
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  double norm = 1.0/(Nx * Ny);
  
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) { 
      // Anti-clockwise plaquette, starting at (x,y). BC handled in read accessor
      plaq += (gauge->read(x,y,0) * gauge->read((x+1),y,1) * conj(gauge->read(x,(y+1),0)) * conj(gauge->read(x,y,1)));
    }
  }
  return plaq * norm;
}

//-----------------------------------------------------------------------------------
