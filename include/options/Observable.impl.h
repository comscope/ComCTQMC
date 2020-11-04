#ifndef INCLUDE_OPTIONS_OBSERVABLE_IMPL_H
#define INCLUDE_OPTIONS_OBSERVABLE_IMPL_H


#include "Observable.impl.h"


namespace opt {

    template<typename Tensor>
    Observable::Observable(Tensor const& tensor) :
    N_(tensor.N()), one_body_(N_*N_, .0), two_body_(N_*N_*N_*N_, .0) {
        for(int fDagg = 0; fDagg < N_; ++fDagg)
        for(int f = 0; f < N_; ++f)
        one_body_[N_*fDagg + f] = tensor.t(fDagg, f);
        
        for(int f1Dagg = 0; f1Dagg < N_; ++f1Dagg)
        for(int f1 = 0; f1 < N_; ++f1)
        for(int f2Dagg = 0; f2Dagg < N_; ++f2Dagg)
        for(int f2 = 0; f2 < N_; ++f2)
        two_body_[N_*N_*N_*f1Dagg +
                  N_*N_*f1 +
                  N_*f2Dagg +
                  f2] = tensor.V(f1Dagg, f1, f2Dagg, f2);
    }


};

#endif //OPTIONS


