#pragma once

#include "schwinger2d_internal.h"
#include "lattice.h"

// strip path from __FILE__
inline constexpr const char* str_end(const char *str) { return *str ? str_end(str + 1) : str; }
inline constexpr bool str_slant(const char *str) { return *str == '/' ? true : (*str ? str_slant(str + 1) : false); }
inline constexpr const char* r_slant(const char* str) { return *str == '/' ? (str + 1) : r_slant(str - 1); }
inline constexpr const char* file_name(const char* str) { return str_slant(str) ? r_slant(str_end(str)) : str; }

//Vector utilities
//---------------------------------------------------------------------------

namespace blas {

  /**
     @brief asserts that two vectors have the same length   
  */
  template <typename T> const void assertVectorLength(const field<T> *x, const field<T> *y,
						      const char *func){
    if(x->size() != y->size()) {
      cout << "Error: vector sizes not equal (" << func << ")" << endl;
      exit(0);
    }
  }
  
  // Zero vector
  template <typename T> void zero(field<T> *x) {
#pragma omp parallel for
    for(int i=0; i<(int)x->size(); i++) x->data[i] = 0.0;
  }
  
  // Copy vector 
  template <typename T> void copy(field<T> *x, const field<T> *y) {
    assertVectorLength(x, y, __func__);
#pragma omp parallel for
    for(int i=0; i<(int)x->size(); i++) x->data[i] = y->data[i];
  }
  
  // Norm squared 
  template <typename T> double norm2(const field<T> *x) {
    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for(int i=0; i<(int)x->size(); i++) sum += (conj(x->data[i]) * x->data[i]).real();
    return sum;
  }
  
  // axpby
  template <typename T, typename A> void axpby(const A a, const field<T> *x, const A b, field<T> *y) {
    assertVectorLength(x,y,__func__);
#pragma omp parallel for
    for(int i=0; i<(int)x->size(); i++) {
      y->data[i] *= b;
      y->data[i] += a*x->data[i];
    }
  }

  // axpy in result
  template <typename T, typename A> void axpy(const A a, const field<T> *x, const field<T> *y, field<T> *z) {
    assertVectorLength(x, y, __func__);
#pragma omp parallel for
    for(int i=0; i<(int)z->size(); i++) {
      z->data[i] = y->data[i] + a*x->data[i];
    }
  }
  
  // axpy in place
  template<typename T, typename A> void axpy(const A a, field<T> *x, field<T> *y) {    
    assertVectorLength(x, y, __func__);
#pragma omp parallel for
    for(int i=0; i<(int)x->size(); i++) {
      y->data[i] += a*x->data[i];
    }
  }
  
  // ax
  template <typename T, typename A> void ax(const A a, field<T> *x) {
#pragma omp parallel for
    for(int i=0; i<(int)x->size(); i++) {
      x->data[i] *= a;
    }
  }

  // Complex inner product
  Complex cDotProd(const field<Complex> *x, const field<Complex> *y);
    
  // Real inner product
  double DotProd(const field<double> *x, const field<double> *y);
  
  // Print the vector elements
  void printVector(const std::vector<Complex> &x);

}
