#ifndef INCLUDE_IO_TENSOR_IMPL_H
#define INCLUDE_IO_TENSOR_IMPL_H

#include "Tensor.h"

namespace io {

template<typename T>
Tensor<T>::Tensor(std::size_t I, std::size_t J, std::size_t K, std::size_t L) : I_(I), J_(J), K_(K), L_(L) {};

template<typename T>
Tensor<T>::Tensor(Tensor&& other) noexcept : I_(other.I_), J_(other.J_),K_(other.K_), L_(other.L_), ijkl_(std::move(other.ijkl_)), data_(std::move(other.data_)), entries_(std::move(other.entries_)) { other.I_ = other.J_ = other.K_ = other.L_ = 0;};

template<typename T>
Tensor<T>& Tensor<T>::operator=(Tensor&& other) { I_ = other.I_; J_ = other.J_; K_ = other.K_; L_ = other.L_; other.I_ = other.J_ = other.K_ = other.L_ = 0; data_ = std::move(other.data_); entries_ = std::move(other.entries_); ijkl_ = std::move(other.ijkl_);  return *this;};

template<typename T>
void Tensor<T>::emplace(int const i, int const j, int const k, int const l, std::string const& entry, T const v){
    
    ijkl_.push_back(std::vector<int>({this->index(i,j,k,l),i,j,k,l}));
    this->emplace(this->index(i,j,k,l),entry,v);
}

template<typename T>
Tensor<T>& Tensor<T>::resize(int I, int J, int K, int L, T value) {
    return *this;
};

template<typename T>
Tensor<T>& Tensor<T>::conj() {
    Tensor temp;
    temp.I_ = J_; temp.J_ = I_; temp.K_ = K_; temp.L_ = L_;
    
    for(int i = 0; i < I_; ++i)
    for(int j = 0; j < J_; ++j)
    for(int k = 0; k < K_; ++k)
    for(int l = 0; l < L_; ++l){
        temp(l, k, j, i) = ut::conj(operator()(i, j, k, l));
    }
    
    return *this = std::move(temp);
};


};

#endif 
