#pragma once

#include "schwinger2d_internal.h"
#include "utils.h"
#include "blas.h"
#include "inverters.h"

//-----------------------------------------------------------------------------------
// 2 Dimensional routines 
//-----------------------------------------------------------------------------------

//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
void measWilsonLoops(field<Complex> *gauge, double plaq, int iter);

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
void measPionCorrelation(field<Complex> *gauge, int iter);

double measTopCharge(field<Complex> *gauge);

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


//-----------------------------------------------------------------------------------
