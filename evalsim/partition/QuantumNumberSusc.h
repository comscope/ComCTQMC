#ifndef EVALSIM_PARTITION_QUANTUMNUMBERSUSC_H
#define EVALSIM_PARTITION_QUANTUMNUMBERSUSC_H


#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/options/Options.h"
#include "../../include/atomic/Generate.h"

namespace evalsim {
    
    namespace partition {
        
        
    jsx::value get_qn_susc(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, jsx::value const& jScalar);
        
    }
}


#endif









