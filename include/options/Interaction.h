#ifndef INCLUDE_OPTIONS_INTERACTION_H
#define INCLUDE_OPTIONS_INTERACTION_H

#include "sphericalharmonics/Basis.h"
#include "sphericalharmonics/SlaterCondon.h"

#include "model/Basis.h"
#include "model/Kanamori.h"

#include "../JsonX.h"
#include "../io/Vector.h"

namespace opt {
    
    struct Interaction {
        Interaction() = delete;
        template<typename Type>
        Interaction(Type const& type);
        
        Interaction(Interaction const& ) = delete;
        Interaction(Interaction&& ) = default;
        Interaction& operator=(Interaction const& ) = delete;
        Interaction& operator=(Interaction&& ) = delete;
        ~Interaction() = default;
        
        double operator()(int m1, int m2, int m3, int m4) const;
        
        int n() const;
        
        std::string approximation() const;
        
    private:
        int const n_;
        std::string const approximation_;
        std::vector<double> tensor_;
    };
    
    
    
    Interaction get_interaction(jsx::value jBasis, jsx::value jTwoBody);
    
};

#include "Interaction.impl.h"

#endif //OPTIONS


