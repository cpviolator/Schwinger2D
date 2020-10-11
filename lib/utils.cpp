#include "utils.h"

template<class T> field<T>::field(param_t p) : p(p)
{
  data.resize(p.Nx * p.Ny * 2);
  for(unsigned int i=0; i<data.size(); i++) data[i] = 0.0;
}

template<class T> field<T>::field(std::vector<T> &data, param_t p) : data(data), p(p) { }    

template<class T>  T field<T>::read(int x, int y, int mu) {
  return data[2*(x%p.Nx + (p.Nx * (y%p.Ny))) + mu];
}
  
template<class T> void field<T>::write(int x, int y, int mu, const T elem) {
    data[2*(x%p.Nx + (p.Nx * (y%p.Ny))) + mu] = elem;
}

template<class T> void field<T>::copy(field<T> *in) {
  blas::copy(data, in->data);
}


template<class T> void field<T>::print() {
  for(int x=0; x<p.Nx; x++) {
    for(int y=0; y<p.Ny; y++) {
      for(int mu=0; mu<2; mu++) {      
	cout << "elem("<<x<<","<<y<<":" << mu << ") = " << *this->read(x,y,mu) << endl;
      }
    }
  }
}

template<class T> field<T>::~field() {
  data.resize(0);
}


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

