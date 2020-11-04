#ifndef INCLUDE_IO_MATRIX_IMPL_H
#define INCLUDE_IO_MATRIX_IMPL_H

#include "Matrix.h"

namespace io {

template <typename T>
Matrix<T>::Matrix(std::size_t I, std::size_t J, T value) : I_(I), J_(J), data_(I_*J_, value) {};

template <typename T>
Matrix<T>::Matrix(Matrix&& other) noexcept : I_(other.I_), J_(other.J_), data_(std::move(other.data_)) { other.I_ = other.J_ = 0;};

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix&& other) { I_ = other.I_; J_ = other.J_; other.I_ = other.J_ = 0; data_ = std::move(other.data_);  return *this;};

template <typename T>
Matrix<T>& Matrix<T>::resize(int I, int J, T value) {
    I_ = I; J_ = J; data_.resize(I_*J_, value);
    return *this;
};

template <typename T>
Matrix<T>& Matrix<T>::conj() {
    Matrix temp; temp.data_.resize(I_*J_);
    temp.I_ = J_; temp.J_ = I_;
    
    for(int i = 0; i < I_; ++i)
    for(int j = 0; j < J_; ++j)
    temp(j, i) = ut::conj(operator()(i, j));
    
    return *this = std::move(temp);
};

template <typename T>
Matrix<T> Matrix<T>::conj() const {
    Matrix temp; temp.data_.resize(I_*J_);
    temp.I_ = J_; temp.J_ = I_;
    
    for(int i = 0; i < I_; ++i)
    for(int j = 0; j < J_; ++j)
    temp(j, i) = ut::conj(operator()(i, j));
    
    return temp;
};

template <typename T>
void Matrix<T>::read(jsx::value const& source) {
    I_ = source(0).int64();
    J_ = source(1).int64();
    data_.read(source(2));
};

template <typename T>
void Matrix<T>::write(jsx::value& dest) const {
    dest = jsx::array_t(3);
    
    dest(0) = static_cast<jsx::int64_t>(I());
    dest(1) = static_cast<jsx::int64_t>(J());
    data_.write(dest(2));
};


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Patch    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



template <typename T>
PrettyMatrix<T>::PrettyMatrix(std::size_t I, std::size_t J, T value) : I_(I), J_(J), data_(I_*J_, value) {};

template <typename T>
PrettyMatrix<T>::PrettyMatrix(PrettyMatrix&& other) noexcept : I_(other.I_), J_(other.J_), data_(std::move(other.data_)) { other.I_ = other.J_ = 0;};

template <typename T>
PrettyMatrix<T>& PrettyMatrix<T>::operator=(PrettyMatrix&& other) { I_ = other.I_; J_ = other.J_; other.I_ = other.J_ = 0; data_ = std::move(other.data_);  return *this;};

template <typename T>
PrettyMatrix<T>& PrettyMatrix<T>::resize(int I, int J, T value) {
    I_ = I; J_ = J; data_.resize(I_*J_, value);
    return *this;
};



};

#endif 
