#ifndef EVALSIM_PARTITION_READDENSITYMATRIX_H
#define EVALSIM_PARTITION_READDENSITYMATRIX_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>


#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"


namespace evalsim {
    
    namespace partition {
        
        namespace meas {

            template<typename Value> // one coud get rid of this fuck when properly implementing meas::Matrix<Value> ....
            jsx::value read_matrix(io::Vector<Value> const& source, int const dim);
            
            template<typename Value>
            jsx::value read_density_matrix(jsx::value const& jParams, jsx::value const& jDensityMatrixMeas);
            
        }
        
    }
    
}

#endif
