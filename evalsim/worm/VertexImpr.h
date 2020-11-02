#ifndef EVALSIM_WORM_VERTEX_IMPR
#define EVALSIM_WORM_VERTEX_IMPR

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Common.h"
#include "include/Vertex.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::vertex_impr::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
    }
    
}

#include "VertexImpr.impl.h"

#endif









