#include "blas.h"

namespace blas {

  /**
     @brief asserts that two vectors have the same length   
  */
  void assertVectorLength(std::vector<Complex> &x, const std::vector<Complex> &y,
			  const char *func){
    if(x.size() != y.size()) {
      cout << "Error: vector sizes not equal (" << __func__ << ")" << endl;
      exit(0);
    }
  }

  // Zero vector
  template<typename T> void zero(std::vector<T> &x) {
    for(int i=0; i<(int)x.size(); i++) x[i] = 0.0;
  }

  // Copy vector 
  template<typename T> void copy(std::vector<T> &x, const std::vector<T> &y) {
    assertVectorLength(x,y,__func__);
    for(int i=0; i<(int)x.size(); i++) x[i] = y[i];
  }

  // Inner product
  template<typename T> T dotProd(const std::vector<T> &x, const std::vector<T> &y) {
    T prod = 0.0;
    assertVectorLength(x,y,__func__);
    for(int i=0; i<(int)x.size(); i++) prod += conj(x[i]) * y[i];
    return prod;
  }
  
  // Norm squared 
  template<typename T> T norm2(std::vector<T> &x) { 
    T norm2 = 0.0;
    for(int i=0; i<(int)x.size(); i++) {
      norm2 += conj(x[i]) * x[i];
    }
    return norm2;
  }

  // Norm 
  double norm(std::vector<Complex> &a) { 
    return sqrt(real(norm2(a)));
  }

  // (c)axpby
  template<typename Ta, typename Tb> void caxpby(const Ta a, std::vector<Complex> &x, const Tb b, std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
    for(int i=0; i<(int)x.size(); i++) {
      y[i] *= b;
      y[i] += a*x[i];
    }
  }

  // (c)axpy in place
  template<typename T> void caxpy(const T a, std::vector<Complex> &x, std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
    for(int i=0; i<(int)x.size(); i++) {
      y[i] += a*x[i];
    }
  }

  // (c)axpy in result
  template<typename T> void caxpy(const T a, std::vector<Complex> &x, std::vector<Complex> &y, std::vector<Complex> &z) {
    assertVectorLength(x,y,__func__);
    for(int i=0; i<(int)x.size(); i++) {
      z[i] = y[i] + a*x[i];
    }
  }

  
  // (c)ax
  template<typename Ta> void cax(const Ta a, std::vector<Complex> &x) {
    for(int i=0; i<(int)x.size(); i++) {
      x[i] *= a;
    }
  }

  // Print the vector elements
  void printVector(const std::vector<Complex> &x){
    for(int i=0; i<(int)x.size(); i++)
      cout << "elem["<<i<<"] = ("<< x[i].real() << "," << x[i].imag() <<")" << endl; 
  }

}
