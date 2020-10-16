#ifndef EVALSIM_PARTITION_SCALAR_H
#define EVALSIM_PARTITION_SCALAR_H


#include "ReadDensityMatrix.h"
#include "ReadHamiltonian.h"


#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/options/Options.h"
#include "../../include/atomic/Generate.h"

namespace evalsim {
    
    namespace partition {
        
        template<typename Value>
        jsx::value get_scalar(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
        
    }
}


#endif









