#ifndef EVALSIM_WORM_VERTEX_IMPRSUM_IMPL
#define EVALSIM_WORM_VERTEX_IMPRSUM_IMPL

#include "VertexImprSum.h"

namespace evalsim {
    
    namespace worm {
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::vertex_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::vertex_imprsum::Worm::name();
        jsx::value jWorm = jParams(name);
        
        return func::vertex::impr_est_evalsim<Value>(jParams, jWorm, jMeasurements, jPartition, jObservables);
    }
    
        
    }
    
}


#endif









