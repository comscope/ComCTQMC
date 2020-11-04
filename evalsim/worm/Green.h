#ifndef EVALSIM_WORM_GREEN
#define EVALSIM_WORM_GREEN

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Green.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"


namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
    
    }
    
}

#include "Green.impl.h"

#endif









