#ifndef EVALSIM_WORM_IMPL_GREEN
#define EVALSIM_WORM_IMPL_GREEN

#include "Green.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::green::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::green::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::green::non_impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
        
    }
    
}


#endif









