#ifndef EVALSIM_PARTITION_OCCUPATIONSUSC_H
#define EVALSIM_PARTITION_OCCUPATIONSUSC_H


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
        io::rmat get_occupation_susc_moments(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, jsx::value& jObservables);
        
        
        template<typename Value>
        jsx::value get_occupation_susc_bulla(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, io::rmat const& moments, io::Matrix<Value> const& occupation, io::Matrix<Value> const& correlation, jsx::value& jObservables);
        
        
        template<typename Value>
        jsx::value get_occupation_susc_direct(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, io::rmat const& moments, io::Matrix<Value> const& occupation, io::Matrix<Value> const& correlation, jsx::value& jObservables);
        
    }
    
}


#endif









