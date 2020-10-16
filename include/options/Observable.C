#include "Observable.h"

namespace opt {


    int Observable::N() const {
        return N_;
    };
    
    double Observable::t(int fDagg, int f) const {
        return one_body_[N_*fDagg + f];
    };
    
    double Observable::V(int f1Dagg, int f1, int f2Dagg, int f2) const {
        return two_body_[N_*N_*N_*f1Dagg +
                         N_*N_*f1 +
                         N_*f2Dagg +
                         f2];
    };
    
    Observable get_observable(jsx::value const& jBasis, std::string const name) {
        if(jBasis("orbitals").is<jsx::string_t>()) {

            if(jBasis("type").string() == "product real" || jBasis("type").string() == "product") {
                if(name == "S2")
                    return Observable(sphericalharmonics::S2(sphericalharmonics::Basis(jBasis)));
            } else if(jBasis("type").string() == "product imag") {
                if(name == "S2")
                    return Observable(sphericalharmonics::S2(sphericalharmonics::Basis(jBasis)));
                //else if(obs.first == "L2")
                //    return Observable(sphericalharmonics::L2(l), transformation);
            } else if(jBasis("type").string() == "coupled") {
                if(name == "J2")
                    return Observable(sphericalharmonics::J2(sphericalharmonics::Basis(jBasis)));
            }
            
            throw std::runtime_error("opt: observable " + name + " not defined for spherical harmonics");

        } else if(jBasis("orbitals").is<jsx::int64_t>()) {

            if(name == "S2")
                return Observable(model::S2(model::Basis(jBasis)));
            
            throw std::runtime_error("opt: observable " + name + " not defined for model basis");
            
        } else
            
            throw std::runtime_error("opt: orbitals not defined");
    }
    
};
