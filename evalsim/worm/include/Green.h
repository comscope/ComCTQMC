#ifndef EVALSIM_WORM_INCLUDE_GREEN
#define EVALSIM_WORM_INCLUDE_GREEN

#include "../../partition/Functions.h"
#include "../../partition/ReadFunctions.h"
#include "../../partition/ReadDensityMatrix.h"
#include "../../partition/ReadHamiltonian.h"

#include "../../../include/JsonX.h"
#include "../../../ctqmc/include/config/Worms.h"

#include "functions/Functions.h"
#include "functions/Measurements.h"
#include "functions/Utilities.h"
#include "Common.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
        
            namespace green {
                
                template<typename Value>
                jsx::value non_impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
            
                template<typename Value>
                jsx::value impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
                
                template<typename Value>
                void compute_green_from_improved(jsx::value const& jParams, iOmega const& iomega, io::Matrix<Value> const& oneBody, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& sigmagreen, std::vector<io::cmat> & green);
                
                template<typename Value>
                void compute_self_from_improved(jsx::value const& jParams, std::vector<io::cmat> const& sigmagreen, std::vector<io::cmat> const& green, std::vector<io::cmat>& self);
                
                template <typename Value>
                std::vector<io::Matrix<Value>> compute_green_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, jsx::value const& jPartition, jsx::value const& jScalar);
                
                template <typename Value>
                std::vector<io::Matrix<Value>> compute_self_moments(jsx::value const& jParams, std::vector<io::Matrix<Value>> const& hybMoments, std::vector<io::Matrix<Value>>& greenMoments);
                
                template <typename Value>
                void add_green_tail(jsx::value const& jParams, iOmega const& iomega, io::Matrix<Value> const& oneBody, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& self, std::vector<io::cmat>& green);
                
                template<typename Value>
                void add_self_tail(jsx::value const& jHybMatrix, iOmega const& iomega, std::vector<io::cmat>& function, std::vector<io::Matrix<Value>> const& moments, std::size_t hybSize);
            }
        }
        
    }
    
}


#endif









