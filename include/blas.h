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
  void assertVectorLength(std::vector<Complex> &x, const std::vector<Complex> &y, const char *func);

  // Zero vector
  template<typename T> void zero(std::vector<T> &x);

  // Copy vector 
  template<typename T> void copy(std::vector<T> &x, const std::vector<T> &y);

  // Inner product
  template<typename T> T dotProd(std::vector<T> &x, const std::vector<T> &y);
  
  // Norm squared 
  double norm2(std::vector<Complex> &x) ;

  // Norm 
  double norm(std::vector<Complex> &a) ;

  // (c)axpby
  template<typename Ta, typename Tb> void caxpby(const Ta a, std::vector<Complex> &x, const Tb b, std::vector<Complex> &y);

  // (c)axpy in place
  template<typename Ta> void caxpy(const Ta a, std::vector<Complex> &x, std::vector<Complex> &y);

  // (c)axpy in result
  template<typename Ta> void caxpy(const Ta a, std::vector<Complex> &x, std::vector<Complex> &y, std::vector<Complex> &z);
  
  // (c)ax
  template<typename Ta> void cax(const Ta a, std::vector<Complex> &x);

  // Print the vector elements
  void printVector(const std::vector<Complex> &x);

}
