#ifndef EVALSIM_WORM_INCLUDE_VERTEX
#define EVALSIM_WORM_INCLUDE_VERTEX

#include "../../partition/Functions.h"
#include "../../partition/ReadFunctions.h"

#include "../../../include/JsonX.h"
#include "../../../ctqmc/include/config/Worms.h"
#include "functions/Functions.h"
#include "functions/Measurements.h"
#include "functions/Utilities.h"
#include "Common.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
        
            namespace vertex {
                
                template<typename Value>
                jsx::value non_impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
            
                template<typename Value>
                jsx::value impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
            
            
                template <typename Value>
                void compute_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::cmat> const& green, OmegaMap const& green_OM, std::vector<io::ctens>& disconnected);
                
                
                template<typename Value>
                void compute_and_subtract_disconnected(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::cmat> const& green, OmegaMap const& green_OM, std::vector<io::ctens>& full_in_connected_out);
                
                
                template<typename Value>
                void compute_connected_from_improved_estimator(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                                               std::vector<io::cmat> const& green,  std::vector<io::cmat> const& self, OmegaMap const& green_OM,
                                                               std::vector<io::ctens> const& improved_estimator, std::vector<io::ctens>& connected);
                
                template<typename Value>
                void enforce_symmetries(jsx::value const& jParams, Frequencies<Value> const& frequencies, std::vector<io::ctens> const& vertex, std::vector<io::ctens>& symm_vertex);
                
                template<typename Value>
                void compute_full_vertex(jsx::value const& jParams, Frequencies<Value> const& frequencies,
                                         std::vector<io::cmat> const& green, OmegaMap const& green_OM,
                                         std::vector<io::ctens> const& susceptibility, std::vector<io::ctens>& full_vertex);
            }
        }
    }
}


#endif









