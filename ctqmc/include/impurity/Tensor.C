
#include "Tensor.h"

namespace imp {
    
    template<typename Value>
    Tensor<Value>::Tensor(jsx::value jTwoBody, int N) :
    N_(N),
    tensor_(jsx::at<io::Vector<Value>>(jTwoBody)),
    non_zero_(N_) {
        if(N_*N_*N_*N_ != tensor_.size())
            throw std::runtime_error("imp::Tensor: interaction tensor has wrong size");
        
        io::Vector<Value> temp(tensor_.size(), .0);
        for(int i = 0; i < N_; ++i)
            for(int j = 0; j < N_; ++j)
                for(int k = 0; k < N_; ++k)
                    for(int l = 0; l < N_; ++l) {
                        auto const& element = ((*this)(i, j, k, l) - (*this)(j, i, k, l) - (*this)(i, j, l, k) + (*this)(j, i, l, k))/4.;
                        
                        if(std::abs(element) > 1.e-14) {
                            temp[N_*N_*N_*i + N_*N_*j + N_*k + l] = element;
                            non_zero_[i].push_back({2*j + 1, 2*k, 2*l});
                        }
                    }
        
        tensor_ = temp;
    };
    
    template struct Tensor<double>;
    template struct Tensor<ut::complex>;
}

