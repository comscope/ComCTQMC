
#include "Tensor.h"

namespace ga {
    

    template<typename Value>
    Tensor<Value>::Tensor(jsx::value jTensor) :
    one_body_(std::move(jsx::at<io::PrettyMatrix<Value>>(jTensor("one body")))),
    two_body_(std::move(jsx::at<io::Vector<Value>>(jTensor("two body")))),
    N_(one_body_.I()) {
        if(one_body_.I() != one_body_.J())
            throw std::runtime_error("ga::Tensor: one body matrix is not square");
        
        if(N_*N_*N_*N_ != static_cast<int>(two_body_.size()))
            throw std::runtime_error("Tensor: one and two body tensor dimensions not compatible");
        
        for(auto& entry : two_body_) entry /= 2.;
    };

    template<typename Value>
    Tensor<Value>::Tensor(Tensor const& tensor, Interaction) :
    one_body_(tensor.N_, tensor.N_),
    two_body_(tensor.two_body_),
    N_(tensor.N_) {
    };
    
    template struct Tensor<std::complex<double>>;
    template struct Tensor<double>;

};






