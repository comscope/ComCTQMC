#ifndef EVALSIM_PARTITION_EVALSIM
#define EVALSIM_PARTITION_EVALSIM

#include <memory>
#include <algorithm>
#include <complex>

#include "ReadFunctions.h"
#include "Functions.h"

#include "Moments.h"
#include "Scalar.h"
#include "Occupation.h"
#include "Probabilities.h"
#include "QuantumNumberSusc.h"
#include "OccupationSusc.h"

#include "../../include/atomic/Generate.h"
#include "../../include/linalg/Operators.h"
#include "../../include/options/Options.h"

#include "../../ctqmc/include/config/Worms.h"

namespace evalsim {
    
    namespace partition {
        
        template<typename Value>
        jsx::value evalsim(jsx::value jParams, jsx::value const& jMeasurements);
        
    }
    
}

    
#endif
    
    
    
    
    
    
    
    
    
