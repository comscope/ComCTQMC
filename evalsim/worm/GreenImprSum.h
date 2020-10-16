#ifndef EVALSIM_WORM_GREEN_IMPRSUM
#define EVALSIM_WORM_GREEN_IMPRSUM

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

#include "../partition/Functions.h"
#include "../partition/ReadFunctions.h"
#include "../partition/ReadDensityMatrix.h"
#include "../partition/ReadHamiltonian.h"

#include "include/Green.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green_imprsum::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
        
    }
    
}

#include "GreenImprSum.impl.h"

#endif









