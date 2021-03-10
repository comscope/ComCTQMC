#ifndef EVALSIM_WORM_INCLUDE_SUSC
#define EVALSIM_WORM_INCLUDE_SUSC

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
            
            namespace susc {
                
                namespace ph {
                
                    template<typename Value>
                    jsx::value non_impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
                
                    template<typename Value>
                    jsx::value impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
                    
                    template <typename Value>
                    void compute_and_subtract_disconnected(jsx::value const& jParams, jsx::value const& jOccupation, BosonFrequencies<Value> const& frequencies, std::vector<io::ctens>& full_in_connected_out);
                
                    template <typename Value>
                    void enforce_symmetries(jsx::value const& jParams, jsx::value const& jWorm, std::vector<io::ctens> const& no_symm, std::vector<io::ctens>& symm);
                    
                    jsx::value qn_susc(jsx::value const& jParams, std::vector<io::ctens> const& susc_tensor);
                
                }
                
                namespace pp {
                
                    template<typename Value>
                    jsx::value non_impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
                
                    template<typename Value>
                    jsx::value impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables);
                    
                    template <typename Value>
                    void compute_and_subtract_disconnected(jsx::value const& jParams, jsx::value const& jOccupation, std::vector<io::ctens>& full_in_connected_out);
                    
                    template <typename Value>
                    void enforce_symmetries(jsx::value const& jParams, std::vector<io::ctens> const& no_symm, std::vector<io::ctens>& symm);
                    
                }
                
            }
            
        }
        
    }
    
}


#endif









