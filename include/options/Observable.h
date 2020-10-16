#ifndef INCLUDE_OPTIONS_OBSERVABLE_H
#define INCLUDE_OPTIONS_OBSERVABLE_H


#include "sphericalharmonics/Observables.h"
#include "model/Observables.h"

#include "../JsonX.h"
#include "../io/Vector.h"


namespace opt {

    struct Observable {
        
        Observable() = delete;
        
        template<typename Tensor>
        Observable(Tensor const& tensor);
        
        Observable(Observable const&) = delete;
        Observable(Observable&& ) = default;
        Observable& operator=(Observable const& ) = delete;
        Observable& operator=(Observable&& ) = default;
        ~Observable() = default;

        int N() const;
        
        double t(int fDagg, int f) const;
        
        double V(int f1Dagg, int f1, int f2Dagg, int f2) const;
        
    private:
        int N_;
        std::vector<double> one_body_;
        std::vector<double> two_body_;
    };
    
    
    Observable get_observable(jsx::value const& jBasis, std::string const name);
    
};

#include "Observable.impl.h"

#endif //OPTIONS


