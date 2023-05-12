#include "measurements.h"

//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
void measWilsonLoops(field<Complex> *gauge, double plaq, int iter)
{
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
  //#pragma omp parallel for
  for(int Xrect=1; Xrect<loop_max; Xrect++) {
      
    //Loop over all Y side sizes of rectangle
    for(int Yrect=1; Yrect<loop_max; Yrect++) {
      
      //Loop over all x,y starting points
      for(x=0; x<Nx; x++)
	for(y=0; y<Ny; y++){
	    
	  Complex w = Complex(1.0,0.0);
	    
	  //Move in +x up to p1.
	  for(int dx=0; dx<Xrect; dx++)     w *= smeared->read((x+dx)%Nx,y,0);
	    
	  //Move in +y up to p2 (p1 constant)
	  int p1 = (x + Xrect)%Nx;
	  for(int dy=0; dy<Yrect; dy++)     w *= smeared->read(p1,(y+dy)%Ny,1);
	  
	  //Move in -x from p1 to (p2 constant)
	  int p2 = (y + Yrect)%Ny;
	  for(int dx=Xrect-1; dx>=0; dx--)  w *= conj(smeared->read((x+dx)%Nx,p2,0));
	  
	  //Move in -y from p2 to y
	  for(int dy=Yrect-1; dy>=0; dy--)  w *= conj(smeared->read(x,(y+dy)%Ny,1));
	  wLoops[Xrect][Yrect] += w*inv_Lsq;
	}
    }
  }

  //Compute string tension
  //#pragma	omp parallel for
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
  snprintf(fname, 100, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d %.16e ", iter, -log(abs(plaq)) );
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
      snprintf(fname, 100, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %.16e %.16e\n", iter, real(wLoops[sizex][sizey]), imag(wLoops[sizex][sizey]));	    
      fclose(fp);
    }
  delete smeared;
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

void measPionCorrelation(field<Complex> *gauge, int iter)
{  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;

  field<Complex> *propUp, *propDn, *source, *Dsource;
  inverterCG *inv = new inverterCG(gauge->p);  
  
  propUp = new field<Complex>(gauge->p);
  propDn = new field<Complex>(gauge->p);
  source = new field<Complex>(gauge->p);
  Dsource = new field<Complex>(gauge->p);
  
  //Up type prop
  blas::zero(source);
  blas::zero(Dsource);
  blas::zero(propUp);
  source->write(0, 0, 0, cUnit);
  
  // up -> (g3Dg3) * up *****
  // (g3Dg3D)^-1 * (g3Dg3) up = D^-1 * up *****
  
  g3psi(Dsource, source);
  g3Dpsi(source, Dsource, gauge);
  inv->solve(propUp, source, gauge);

  //Down type prop
  blas::zero(source);
  blas::zero(Dsource);
  blas::zero(propDn);
  source->write(0, 0, 1, cUnit);
  
  // dn -> (g3Dg3) * dn   
  g3psi(Dsource, source);
  g3Dpsi(source, Dsource, gauge);

  // (g3Dg3D)^-1 * (g3Dg3) dn = D^-1 * dn 
  inv->solve(propDn, source, gauge);
  
  double pion_corr[Ny];    
  //Let y be the 'time' dimension
  double corr = 0.0, tmp = 0.0;
  for(int y=0; y<Ny; y++) {
    //initialise
    pion_corr[y] = 0.0;
    //Loop over space and spin, fold propagator
    corr = 0.0;
    for(int x=0; x<Nx; x++) {
      tmp = abs((conj(propDn->read(x,y,0)) * propDn->read(x,y,0) +
      		 conj(propDn->read(x,y,1)) * propDn->read(x,y,1) +
		 conj(propUp->read(x,y,0)) * propUp->read(x,y,0) +
      		 conj(propUp->read(x,y,1)) * propUp->read(x,y,1)));
      
      corr += tmp;
    }
    //NOTE: CHANGE TO UNFOLDED
    //Compute folded propagator
    if ( y < ((Ny/2)+1) ) pion_corr[y] += corr;
    else {
      pion_corr[Ny-y] += corr;
      pion_corr[Ny-y] /= 2.0;
    }
  }
  
  char fname[256];
  string name;
  FILE *fp;
  
  //Full pion correlation
  name = "data/pion/pion";
  constructName(name, gauge->p);
  name += ".dat";  
  snprintf(fname, 100, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter);
  for(int t=0; t<Ny/2+1; t++)
    fprintf(fp, "%.16e ", pion_corr[t]);
  fprintf(fp, "\n");
  fclose(fp);

  delete propUp;
  delete propDn;
  delete source;
  delete Dsource;
  delete inv;
}

double measTopCharge(field<Complex> *gauge){
  
  field<Complex> *smeared = new field<Complex>(gauge->p);
  smearLink(smeared, gauge);

  double top = 0.0;
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  //#pragma	omp parallel for reduction(+:top)
  for(int x=0; x<Nx; x++) {
    int xp1 = (x+1)%Nx;
    for(int y=0; y<Ny; y++){
      int yp1 = (y+1)%Ny;
      Complex w = smeared->read(x,y,0) * smeared->read(xp1,y,1) * conj(smeared->read(x,yp1,0)) * conj(smeared->read(x,y,1));
      top += arg(w)/TWO_PI;  // -pi < arg(w) < pi  Geometric value is an integer.
    }
  }
  delete smeared;
  //return top/TWO_PI;
  return top;
}

//-----------------------------------------------------------------------------------

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

*/


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

