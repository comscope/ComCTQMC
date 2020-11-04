#ifndef EVALSIM_WORM_HEDIN_IMPRSUM
#define EVALSIM_WORM_HEDIN_IMPRSUM

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "../partition/Functions.h"
#include "../partition/ReadFunctions.h"

#include "include/Common.h"
#include "include/Hedin.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_ph_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
    
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_pp_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
    }
    
}


#include "HedinImprSum.impl.h"


#endif









