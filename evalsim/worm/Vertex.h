#ifndef EVALSIM_WORM_VERTEX
#define EVALSIM_WORM_VERTEX

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
        jsx::value evalsim(ut::wrap<cfg::vertex::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
    }
    
}

#include "Vertex.impl.h"

#endif









