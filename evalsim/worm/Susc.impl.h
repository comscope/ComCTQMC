#ifndef EVALSIM_WORM_SUSC_IMPL
#define EVALSIM_WORM_SUSC_IMPL

#include "Susc.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::susc_ph::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::susc_ph::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::susc::ph::non_impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::susc_pp::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::susc_pp::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::susc::pp::non_impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
        
        
    }
    
}


#endif









