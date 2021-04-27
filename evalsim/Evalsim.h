#ifndef EVALSIM_EVALSIM_H
#define EVALSIM_EVALSIM_H


#include "partition/Evalsim.h"
#include "worm/Evalsim.h"

#include "../ctqmc/include/config/Worms.h"
#include "../ctqmc/include/Params.h"
#include "../include/JsonX.h"
#include "../include/measurements/Measurements.h"
#include "../include/measurements/Error.h"
#include "../include/io/Vector.h"
#include "../include/io/Tag.h"
#include "../include/mpi/Utilities.h"
#include "../include/parameters/Initialize.h"

namespace evalsim {
    

//Main interfaces
    void evalsim_driver(const char* case_name);

//---------------------------------------------------------------------------------------------------------------------
    
    template <typename Value>
    void complete_params(jsx::value & jParams);

    jsx::value get_observables(jsx::value & jParams, std::string const name);

    namespace worm {
        
        // catch worm-evalsims that are not yet implemented
        template<typename Value, typename W>
        jsx::value evalsim(ut::wrap<W>, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
            throw std::runtime_error("evalsim for " + W::name() + " worm not implemented");
        }
        
    }
    
    template<typename Value>
    bool add_dynamic(jsx::value& jStatic, jsx::value const& jDynamic);

    void add_dynamics(jsx::value const& jParams, jsx::value& jMeasurements, std::string const worm, std::string const meas);

    template<typename Value>
    struct worm_clean_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, jsx::value& jMeasurements) const {
            if( W::name() != cfg::partition::Worm::name() && jParams.is(W::name()) ) {
                jsx::value temp = std::move(jMeasurements(W::name())("static"));
                
                jMeasurements(W::name()) = std::move(temp);
            }
        }
    };
    
    template<typename Value>
    struct worm_evalsim_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value& jObservables) const {
            if( W::name() != cfg::partition::Worm::name() && jParams.is(W::name()) ) {
                mpi::cout << std::endl << "Begin evaluating " + W::name() + " worm measurements" << std::endl;
                
                jObservables[W::name()] = worm::evalsim<Value>(w, jParams, jMeasurements(W::name()), jMeasurements(cfg::partition::Worm::name()), jObservables);
                
                mpi::cout << "End evaluating " + W::name() + " worm measurements" << std::endl;
            }
        }
        void operator()(ut::wrap<cfg::partition::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jObservables) const {
        };
    };
    
    
    template<typename Value>
    jsx::value evalsim(jsx::value & jParams, jsx::value& jMeasurements);
        
    
}

#endif












