#ifndef CTQMC_INCLUDE_IMPURITY_TENSOR_H
#define CTQMC_INCLUDE_IMPURITY_TENSOR_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../Utilities.h"
#include "../../../include/JsonX.h"
#include "../../../include/io/Vector.h"


namespace imp {
    
    struct TensorIndices {
        TensorIndices() = delete;
        TensorIndices(int flavor1, int flavor2, int flavor3) : flavor1_(flavor1), flavor2_(flavor2), flavor3_(flavor3) {};
        
        int flavor1() const { return flavor1_;};
        int flavor2() const { return flavor2_;};
        int flavor3() const { return flavor3_;};
        
    private:
        int flavor1_, flavor2_, flavor3_;
    };
    
    
    // all this needs to be modified in the presence of superconductivity
    
    
    template<typename Value>
    struct Tensor {
        Tensor() = delete;
        Tensor(jsx::value jTwoBody, int N);
        
        Value operator()(int i, int j, int k, int l) const {
            return tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l];
        };
        
        std::vector<TensorIndices> const& nonzero(int i) const{
            return non_zero_[i];
        };
        
        //std::vector<TensorIndices> const& nonzerodagg(int i) const{
        //    return non_zero_[i];
        //};
        
    private:
        int const N_;
        io::Vector<Value> tensor_;
        std::vector<std::vector<TensorIndices>> non_zero_;
        
        Value& operator()(int i, int j, int k, int l) {
            return tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l];
        };
        
    };
    
}

#include "Tensor.impl.h"

#endif
