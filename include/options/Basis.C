#include "Basis.h"

namespace opt {
    
    template <typename Value>
    int Basis<Value>::n() const {
        return n_;
    };

    template <typename Value>
    int Basis<Value>::N() const {
        return transformation_.I();
    };

    template <typename Value>
    std::complex<double> Basis<Value>::operator()(int f, int m, int s) const {
        return u_[f*2*n_ + 2*m + s];
    };

    template <typename Value>
    std::map<std::string, io::rvec> const& Basis<Value>::qns() const {
        return qns_;
    };

    template<typename Value>
    Basis<Value> get_basis(jsx::value jBasis) {
        jsx::value jTransformation = jBasis.is("transformation") ? jBasis("transformation") : jsx::empty_t();
        
        if(jBasis("orbitals").is<jsx::string_t>()) {
            
            if(jBasis("type").string() == "product real" || jBasis("type").string() == "product")
                return Basis<Value>(sphericalharmonics::Real(jBasis), jTransformation);
            else if(jBasis("type").string() == "product imag")
                return Basis<Value>(sphericalharmonics::Imag(jBasis), jTransformation);
            else if(jBasis("type").string() == "coupled")
                return Basis<Value>(sphericalharmonics::Coupled(jBasis), jTransformation);
            else
                throw std::runtime_error("opt: type " + jBasis("type").string() + " not defined for spherical harmonics basis");
            
        } else if(jBasis("orbitals").is<jsx::int64_t>()) {
            
            return Basis<Value>(model::Basis(jBasis), jTransformation);
            
        } else
            
            throw std::runtime_error("opt: orbitals not defined");
    }
    
    template Basis<double> get_basis<double>(jsx::value jBasis);
    template Basis<std::complex<double>> get_basis<std::complex<double>>(jsx::value jBasis);
    
    template struct Basis<double>;
    template struct Basis<std::complex<double>>;
    
};
