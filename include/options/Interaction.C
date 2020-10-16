#include "Interaction.h"

namespace opt {
    
    double Interaction::operator()(int m1, int m2, int m3, int m4) const {
        return tensor_[n_*n_*n_*m1 +
                       n_*n_*m2 +
                       n_*m3 +
                       m4];
    };
    
    int Interaction::n() const {
        return n_;
    };
    
    std::string Interaction::approximation() const {
        return approximation_;
    };
    
    
    Interaction get_interaction(jsx::value jBasis, jsx::value jTwoBody) {
        if(jBasis("orbitals").is<jsx::string_t>()) {
            
            sphericalharmonics::Basis basis(jBasis);
            
            if(jTwoBody("parametrisation").string() == "slater-condon")
                return Interaction(sphericalharmonics::SlaterCondon(basis, jTwoBody));
            else
                throw std::runtime_error("opt: parametrisation" + jTwoBody("parametrisation").string() + " not defined for spherical harmonics basis");
            
        } else if(jBasis("orbitals").is<jsx::int64_t>()) {
            
            model::Basis basis(jBasis);
            
            if(jTwoBody("parametrisation").string() == "kanamori")
                return Interaction(model::Kanamori(basis, jTwoBody));
            else
                throw std::runtime_error("opt: parametrisation" + jTwoBody("parametrisation").string() + " not defined for model basis");
            
        } else
            throw std::runtime_error("opt: orbitals not defined");
    }
    
    
};

