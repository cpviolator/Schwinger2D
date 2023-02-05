#pragma once

template<typename T> class field {
  
public:
  
  ~field(){
    data.resize(0);
  }
  
  std::vector<T> data;
  Param p;

  int tag;
  
  field(Param p) : p(p)
  {
    data.resize(p.Nx * p.Ny * 2);
    for(unsigned int i=0; i<data.size(); i++) data[i] = 0.0;
  }
  
  T read(int x, int y, int mu) const {
    return data[2*(x + p.Nx * y) + mu];
  } 
  
  void write(int x, int y, int mu, const T elem) {
    data[2*(x + p.Nx * y) + mu] = elem;
  }
  
  unsigned int const size() const { return data.size(); }
  
  T elem(int i) { return data[i]; }
  
  void print() {
    for(int x=0; x<p.Nx; x++) {
      for(int y=0; y<p.Ny; y++) {
	for(int mu=0; mu<2; mu++) {      
	  cout << "elem("<<x<<","<<y<<":" << mu << ") = " << data[2*(x + p.Nx * y) + mu] << endl;
	}
      }
    }
  }
  
  void print(int n) {
    cout << "elem " << n << " = " << data[n] << endl;
  }
  
};
