#ifndef EVALSIM_WORM_VERTEX_IMPRSUM
#define EVALSIM_WORM_VERTEX_IMPRSUM


#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "../partition/Functions.h"
#include "../partition/ReadFunctions.h"

#include "include/Common.h"
#include "include/Vertex.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::vertex_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
    }
    
}

#include "VertexImprSum.impl.h"

#endif









