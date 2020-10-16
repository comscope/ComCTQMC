#ifndef INCLUDE_OPTIONS_BASIS_H
#define INCLUDE_OPTIONS_BASIS_H

#include "Transformation.h"
#include "sphericalharmonics/Basis.h"
#include "model/Basis.h"

#include "../JsonX.h"
#include "../io/Vector.h"


namespace opt {
    
    template<typename Value>
    struct Basis {
        Basis() = delete;
        
        template<typename Type>
        Basis(Type const& type, jsx::value const& jTransformation);
        
        Basis(Basis const&) = delete;
        Basis(Basis&& ) = default;
        Basis& operator=(Basis const& ) = delete;
        Basis& operator=(Basis&& ) = default;
        ~Basis() = default;
        
        int n() const;
        
        int N() const;
        
        std::complex<double> operator()(int f, int m, int s) const;
        
        std::map<std::string, io::rvec> const& qns() const;
        
    private:
        int n_;
        Transformation<Value> const transformation_;
        std::vector<std::complex<double>> u_;
        std::map<std::string, io::rvec> qns_;
    };
    
    
    template<typename Value>
    Basis<Value> get_basis(jsx::value jBasis);
    
};

#include "Basis.impl.h"

#endif //OPTIONS


