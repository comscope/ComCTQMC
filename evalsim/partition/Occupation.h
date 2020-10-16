#ifndef EVALSIM_PARTITION_OCCUPATION_H
#define EVALSIM_PARTITION_OCCUPATION_H


#include "ReadDensityMatrix.h"
#include "Functions.h"

#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/mpi/Utilities.h"

namespace evalsim {
    
    namespace partition {
        
        
        template<typename Value>
        jsx::value get_occupation(jsx::value const& jParams, jsx::value const& jMeasurements, io::Matrix<Value>& occupation, io::Matrix<Value>& correlation);
        
    }
}


#endif









