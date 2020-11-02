#ifndef EVALSIM_WORM_VERTEX_IMPR_IMPL
#define EVALSIM_WORM_VERTEX_IMPR_IMPL

#include "VertexImpr.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::vertex_impr::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::vertex_impr::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::vertex::impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
        }
    
        
    }
    
}


#endif









