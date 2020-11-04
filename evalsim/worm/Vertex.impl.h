#ifndef EVALSIM_WORM_VERTEX_IMPL
#define EVALSIM_WORM_VERTEX_IMPL

#include "Vertex.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::vertex::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::vertex::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::vertex::non_impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
        
    }
    
}


#endif









