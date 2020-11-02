#ifndef EVALSIM_WORM_HEDIN_IMPL
#define EVALSIM_WORM_HEDIN_IMPL

#include "Hedin.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::hedin_ph::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::hedin_ph::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::hedin::ph::non_impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::hedin_pp::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::hedin_pp::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::hedin::pp::non_impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
        
        
    }
    
}


#endif









