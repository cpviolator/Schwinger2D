#include "utils.h"

void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  XSize = "<< p.Nx << endl;
  cout << "          YSize = "<< p.Ny << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Dynamic = " << (p.dynamic == true ? "True" : "False") << endl;
  if (p.dynamic == true) cout << "          Mass = " << p.m << endl;
  cout << "HMC:      Therm Sweeps: (" << p.therm << " accept) (" << p.therm << " accept/reject)" << endl; 
  cout << "          Data Points = " << p.iter_hmc << endl;
  cout << "          Time Step = " << p.tau/p.n_step << endl;
  cout << "          Trajectory Steps " << p.n_step << endl;
  cout << "          Trajectory Length = " << p.tau << endl;
  cout << "Smearing: APE iter = " << p.smear_iter << endl;
  cout << "          APE alpha = " << p.alpha << endl;
#ifdef USE_VOATOL
  cout << "VOATOL:   nkv = " << p.n_kr << endl;
  cout << "          nev = " << p.n_ev << endl;
  cout << "          tol = " << p.eig_tol << endl;
  cout << "          maxiter = " << p.eig_max_restart << endl;
#endif  
}

void constructName(string &name, param_t p) {
  name += "_LX" + to_string(p.Nx) + "_LY" + to_string(p.Ny) + "_B" + to_string(p.beta);
  if(p.dynamic == true) name += "_M"+ to_string(p.m);
  name += "_tau" + to_string(p.tau) + "_nHMCstep" + to_string(p.n_step);
}

void writeGauge(field<Complex> *gauge, string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  //outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge).real() << endl;

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;  
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	outPutFile << setprecision(12) <<  setw(20) << arg(gauge->read(x,y,mu)) << endl;
  
  outPutFile.close();
  return;  
}

void readGauge(field<Complex> *gauge, string name)
{

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  //Header check
  //getline(inPutFile, val);
  //double plaq_real_header = stod(val);
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	getline(inPutFile, val);
	gauge->write(x, y, mu, polar(1.0, stod(val)));
      }
    }
  }
  
  //double plaq_real_measured = measPlaq(gauge).real();
  //double err = fabs(1.0 - plaq_real_header/plaq_real_measured);
  //if(abs(err) > 1e-12) {
  //cout << "Gauge read fail!" << endl;
  //cout << setprecision(16) << setw(20) << "Plaqette on file  = " << plaq_real_header << endl;
  //cout << setprecision(16) << setw(20) << "Plaqette measured = " << plaq_real_measured << endl;   
  //exit(0);
  //}    
  return;
}

void gaussStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	gauge->write(x, y, mu, polar(1.0,drand48()));  
}  

void coldStart(field<Complex> *gauge) {

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int x=0; x<Nx; x++) 
    for(int y=0; y<Ny; y++)
      for(int mu=0; mu<2; mu++)
	gauge->write(x, y, mu, Complex(1.0,0.0));    
}  

// Normalized gaussian exp(-phi*phi/2) and  <phi|phi> = 1
void gaussReal(field<double> *field) {
  
  double r, theta, sum;
  int Nx = field->p.Nx;
  int Ny = field->p.Ny;
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	
	r = sqrt(-2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field->write(x,y,mu, r*cos(theta));
      }
    }
  }  
}

//normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
void gaussComplex(field<Complex> *field) {

  double r1, theta1, r2, theta2, sum;
  double inv_sqrt2 = 1.0/sqrt(2);

  int Nx = field->p.Nx;
  int Ny = field->p.Ny;
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int mu=0; mu<2; mu++) {
	
	r1 = sqrt(-2.0*log(drand48()));
	theta1 = TWO_PI*(drand48());
	r2 = sqrt(-2.0*log(drand48()));
	theta2 = TWO_PI*(drand48());
	
	field->write(x,y,mu, Complex(r1*cos(theta1),r2*sin(theta2))*inv_sqrt2);
      }
    }
  }
}

//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(field<Complex> *smeared, field<Complex> *gauge){

  double alpha = gauge->p.alpha;
  int iter = gauge->p.smear_iter;
  Complex tmp = 0;
  int xp1, xm1, yp1, ym1;
  
  field<Complex> *smeared_tmp = new field<Complex>(gauge->p);
  smeared->copy(gauge);
  smeared_tmp->copy(smeared);
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  for(int i=0; i<iter; i++) {    
    for(int x=0; x<Nx; x++) {
      for(int y=0; y<Ny; y++) {

	xp1 = (x+1)%Nx;
	xm1 = (x-1+Nx)%Nx;
	yp1 = (y+1)%Ny;
	ym1 = (y-1+Ny)%Ny;
		
	//|->-|   |   |
	//^   v + v   ^
	//|   |   |->-|
	tmp = alpha * (smeared->read(x,y,1) * smeared->read(x,yp1,0) * conj(smeared->read(xp1,y,1)));
	
	tmp += alpha * (conj(smeared->read(x,ym1,1)) * smeared->read(x,ym1,0) * smeared->read(xp1,ym1,1));
			
	smeared_tmp->write(x,y,0, tmp);
				
	//|->-    -<-|
	//^    +     ^
	//|-<-    ->-|
	tmp = alpha * (smeared->read(x,y,0) * smeared->read(xp1,y,1) * conj(smeared->read(x,yp1,0)));
	
	tmp += alpha * (conj(smeared->read(xm1,y,0)) * smeared->read(xm1,y,1) * smeared->read(xm1, yp1,01));
	
	smeared_tmp->write(x,y,1, tmp);	
      }
    }
    
    //Project back to U(1)
    for(int x=0; x<Nx; x++)
      for(int y=0; y<Ny; y++)
	for(int mu=0; mu<2; mu++)
	  smeared->write(x,y,mu, polar(1.0,arg(smeared_tmp->read(x,y,mu))));
  }
}
