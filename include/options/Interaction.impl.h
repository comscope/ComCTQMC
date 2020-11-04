#ifndef INCLUDE_OPTIONS_INTERACTION_IMPL_H
#define INCLUDE_OPTIONS_INTERACTION_IMPL_H

#include "Interaction.h"

namespace opt {
    
    template<typename Type>
    Interaction::Interaction(Type const& type) :
    n_(type.n()),
    approximation_(type.approximation()),
    tensor_(n_*n_*n_*n_) {
        for(int m1 = 0; m1 < n_; ++m1)
        for(int m2 = 0; m2 < n_; ++m2)
        for(int m3 = 0; m3 < n_; ++m3)
        for(int m4 = 0; m4 < n_; ++m4)
        tensor_[n_*n_*n_*m1 +
                n_*n_*m2 +
                n_*m3 +
                m4] = type(m1, m2, m3, m4);
    }
        
};

#endif //OPTIONS


