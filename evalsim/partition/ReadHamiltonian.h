#ifndef EVALSIM_PARTITION_READHAMILTONIAN_H
#define EVALSIM_PARTITION_READHAMILTONIAN_H

#include <complex>
#include <vector>
#include <map>


#include "../../include/linalg/Operators.h"
#include "../../include/atomic/Generate.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"

namespace evalsim {
    
    namespace partition {
            
        template<typename Value>
        jsx::value get_hamiltonian(jsx::value const& jParams);
            
            
        template<typename Value>
        jsx::value get_effective_hamiltonian(jsx::value const& jParams);
            
    }
}


#endif









