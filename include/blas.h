#pragma once

#include "schwinger2d_internal.h"

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
  void assertVectorLength(const std::vector<Complex> &x, const std::vector<Complex> &y, const char *func);

  void assertVectorLength(const std::vector<double> &x, const std::vector<double> &y, const char *func);

  // Zero vector
  void zero(std::vector<Complex> &x);
  
  // Zero vector
  void zero(std::vector<double> &x);

  // Copy vector 
  void copy(std::vector<Complex> &x, const std::vector<Complex> &y);

  // Copy vector 
  void copy(std::vector<double> &x, const std::vector<double> &y);
  
  // Inner product
  Complex cDotProd(const std::vector<Complex> &x, const std::vector<Complex> &y);
  
  // Norm squared 
  double norm2(std::vector<Complex> &x);
  
  // Norm squared 
  double norm2(Complex *x, int size);
  
  // Norm 
  double norm(std::vector<Complex> &a);

  // caxpby
  void caxpby(const Complex a, const std::vector<Complex> &x, const Complex b, std::vector<Complex> &y);

  // axpby
  void axpby(const double a, const std::vector<Complex> &x, const double b, std::vector<Complex> &y);

  // caxpy in place
  void caxpy(const Complex a, const std::vector<Complex> &x, std::vector<Complex> &y);

  // caxpy in result
  void caxpy(const Complex a, const std::vector<Complex> &x, const std::vector<Complex> &y, std::vector<Complex> &z);

  // axpy in place
  void axpy(const double a, const std::vector<Complex> &x, std::vector<Complex> &y);

  // axpy in result
  void axpy(const double a, const std::vector<Complex> &x, const std::vector<Complex> &y, std::vector<Complex> &z);

  // axpy in place
  void axpy(const double a, const std::vector<double> &x, std::vector<double> &y);
  
  // axpy in result
  void axpy(const double a, const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &z);

  
  // cax
  void cax(const Complex a, std::vector<Complex> &x);
  
  // ax
  void ax(const double a, std::vector<Complex> &x);

  // cax
  void cax(const Complex a, std::vector<double> &x);
  
  // ax
  void ax(const double a, std::vector<double> &x);

  
  // Print the vector elements
  void printVector(const std::vector<Complex> &x);

}
