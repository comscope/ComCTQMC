#ifndef EVALSIM_PARTITION_MOMENTS_H
#define EVALSIM_PARTITION_MOMENTS_H


#include "ReadDensityMatrix.h"
#include "ReadHamiltonian.h"
#include "Functions.h"

#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"

namespace evalsim {
    
    namespace partition {
        
        
        template<typename Value>
        std::vector<io::Matrix<Value>> get_green_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, jsx::value const& jMeasurements, jsx::value const& jScalar);
        
        
        template<typename Value>
        std::vector<io::Matrix<Value>> get_self_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, std::vector<io::Matrix<Value>> const& greenMoments);
        
    }
}


#endif









