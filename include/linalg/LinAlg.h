#ifndef INCLUDE_LINALG_LINALG_H
#define INCLUDE_LINALG_LINALG_H

#include <vector>
#include <complex>
#include <cassert>
#include <algorithm>

#include "../io/Matrix.h"
#include "../BlasLapack.h"

namespace linalg {
    
template<typename Value>
void mult(char transA, char transB, Value alpha, io::Matrix<Value> const& matrixA, io::Matrix<Value> const& matrixB, Value beta, io::Matrix<Value>& matrixC);
    
template<typename Value>
Value trace(io::Matrix<Value> const& matrix);

template<typename Value>
Value trace(io::Matrix<Value> const& matrixA, io::Matrix<Value> const& matrixB);
    
void eig(char jobz, char uplo, io::rmat& matrix, io::rvec& eigen_values);
void eig(char jobz, char uplo, io::cmat& matrix, io::rvec& eigen_values);

template<typename Value>
io::Matrix<Value> inv(io::Matrix<Value> const& matrix);
    
double spectral_norm(io::rmat const& matrix);
double spectral_norm(io::cmat const& matrix);

template <typename T>
typename std::enable_if<!std::is_integral<T>::value, std::vector<T>>::type linspace(std::size_t const size, T const begin, T const end, bool endpoint = true);
 
}

#include "LinAlg.impl.h"

#endif
