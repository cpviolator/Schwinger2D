#include "blas.h"

namespace blas {

  /**
     @brief asserts that two vectors have the same length   
  */
  void assertVectorLength(const std::vector<Complex> &x, const std::vector<Complex> &y,
			  const char *func){
    if(x.size() != y.size()) {
      cout << "Error: vector sizes not equal (" << __func__ << ")" << endl;
      exit(0);
    }
  }

  void assertVectorLength(const std::vector<double> &x, const std::vector<double> &y,
			  const char *func){
    if(x.size() != y.size()) {
      cout << "Error: vector sizes not equal (" << __func__ << ")" << endl;
      exit(0);
    }
  }

  
  // Zero vector
  void zero(std::vector<Complex> &x) {
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) x[i] = 0.0;
  }

  // Zero vector
  void zero(std::vector<double> &x) {
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) x[i] = 0.0;
  }

  // Copy vector 
  void copy(std::vector<Complex> &x, const std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) x[i] = y[i];
  }

  // Copy vector 
  void copy(std::vector<double> &x, const std::vector<double> &y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) x[i] = y[i];
  }

  // Inner product
  Complex cDotProd(const std::vector<Complex> &x, const std::vector<Complex> &y) {
    Complex prod = 0.0;
    assertVectorLength(x,y,__func__);
#pragma omp parallel for reduction(+:prod)
    for(int i=0; i<(int)x.size(); i++) prod += conj(x[i]) * y[i];
    return prod;
  }
  
  // Norm squared
  double norm2(std::vector<Complex> &x) { 
    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for(int i=0; i<(int)x.size(); i++) sum += (conj(x[i]) * x[i]).real();
    return sum;
  }
  
  // Norm 
  double norm(std::vector<Complex> &a) { 
    return sqrt(real(blas::norm2(a)));
  }
  
  // caxpby
  void caxpby(const Complex a, const std::vector<Complex> &x, const Complex b, std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      y[i] *= b;
      y[i] += a*x[i];
    }
  }

  // axpby
  void axpby(const double a, const std::vector<Complex> &x, const double b, std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      y[i] *= b;
      y[i] += a*x[i];
    }
  }
  
  // caxpy in place
  void caxpy(const Complex a, const std::vector<Complex> &x, std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      y[i] += a*x[i];
    }
  }
  
  // caxpy in result
  void caxpy(const Complex a, const std::vector<Complex> &x, const std::vector<Complex> &y, std::vector<Complex> &z) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      z[i] = y[i] + a*x[i];
    }
  }

  // axpy in place
  void axpy(const double a, const std::vector<Complex> &x, std::vector<Complex> &y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      y[i] += a*x[i];
    }
  }

  // axpy in result
  void axpy(const double a, const std::vector<Complex> &x, const std::vector<Complex> &y, std::vector<Complex> &z) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      z[i] = y[i] + a*x[i];
    }
  }

  // (c)ax
  void cax(const Complex a, std::vector<Complex> &x) {
#pragma omp parallel for
    for(int i=0; i<(int)x.size(); i++) {
      x[i] *= a;
    }
  }

  void ax(const double a, std::vector<Complex> &x) {
#pragma omp parallel for
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
