#ifndef EVALSIM_WORM_SUSC
#define EVALSIM_WORM_SUSC


#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "include/Susc.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::susc_ph::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::susc_pp::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
        
    }
    
}

#include "Susc.impl.h"

#endif









