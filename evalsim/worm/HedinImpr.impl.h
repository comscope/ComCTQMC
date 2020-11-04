#ifndef EVALSIM_WORM_HEDIN_IMPR_IMPL
#define EVALSIM_WORM_HEDIN_IMPR_IMPL

#include "HedinImpr.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::hedin_ph_impr::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::hedin_ph_impr::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::hedin::ph::impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::hedin_pp_impr::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::hedin_pp_impr::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::hedin::pp::impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
        
        
    }
    
}




#endif









