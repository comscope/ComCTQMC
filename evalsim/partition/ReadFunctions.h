#ifndef EVALSIM_PARTITION_READFUNCTIONS_H
#define EVALSIM_PARTITION_READFUNCTIONS_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "Bessel.h"
#include "Functions.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"


namespace evalsim {
    
    namespace partition {
        
        namespace meas {
            
            namespace impl {
                
                template<typename Value>
                std::pair<io::cvec, io::cvec> read_matsubara_basis(io::Vector<Value> const& data, jsx::value const& jParams, jsx::value const& jPartition, int hybSize);
                
                
                template<typename Value>
                std::pair<io::cvec, io::cvec> read_legendre_basis(io::Vector<Value> const& data, jsx::value const& jParams, jsx::value const& jPartition, int hybSize);
                
                
                template<typename Value>
                std::pair<io::cvec, io::cvec> read_function(jsx::value const& jFunction, jsx::value const& jParams, jsx::value const& jPartition, int hybSize);
                
                
                template<typename Value>
                std::pair<std::vector<io::cmat>, std::vector<io::cmat>> read_functions(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize);
                
                
                std::vector<io::cmat> symmetrize_functions(std::vector<io::cmat> const& functions, std::vector<io::cmat> const& functionsConj, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize);
                    
            }
            
            template<typename Value>
            std::vector<io::cmat> read_functions(jsx::value const& jFunctions, jsx::value const& jParams, jsx::value const& jPartition, int hybSize);
            
            template<typename Value>
            std::pair<std::vector<io::cmat>, std::vector<io::cmat>> read_conj_functions(jsx::value const& jFunctions, jsx::value const& jFunctionsConj, jsx::value const& jParams, jsx::value const& jPartition, int hybSize);
            
        }
        
    }
    
}

#endif
