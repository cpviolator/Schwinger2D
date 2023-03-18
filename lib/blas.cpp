#include "blas.h"

namespace blas {

  // Inner product
  Complex cDotProd(const field<Complex> *x, const field<Complex> *y) {
    double re_prod = 0.0;
    double im_prod = 0.0;
    assertVectorLength(x,y,__func__);
#pragma omp parallel for reduction(+:re_prod, im_prod)
    for(int i=0; i<(int)x->size(); i++) {
      Complex prod = conj(x->data[i]) * y->data[i];
      re_prod += prod.real();
      im_prod += prod.imag();
    }
    return Complex(re_prod, im_prod);
  }
    
      // Real Inner product
  double DotProd(const field<double> *x, const field<double> *y) {
    double re_prod = 0.0;
    assertVectorLength(x,y,__func__);
#pragma omp parallel
    for(int i=0; i<(int)x->size(); i++) {
      double prod = x->data[i] * y->data[i];
      re_prod += prod;
    }
    return re_prod;
  }
  
  // Print the vector elements
  void printVector(const std::vector<Complex> &x){
    for(int i=0; i<(int)x.size(); i++)
      cout << "elem["<<i<<"] = ("<< x[i].real() << "," << x[i].imag() <<")" << endl; 
  }  
}
