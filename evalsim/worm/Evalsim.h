#ifndef EVALSIM_WORM_EVALSIM
#define EVALSIM_WORM_EVALSIM


#include "Eval2Ops.h"
#include "Eval4Ops.h"
#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"


namespace evalsim {
    
    namespace worm {
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_2_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }

        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green_impr::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_2_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::green_imprsum::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_2_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::susc_ph::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::susc_pp::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }

        
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_ph::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_ph_impr::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_ph_imprsum::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_pp::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_pp_impr::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::hedin_pp_imprsum::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::vertex::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::vertex_impr::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
        template<typename Value>
        jsx::value evalsim(ut::wrap<cfg::vertex_imprsum::Worm> w, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            return eval_4_ops<Value>(w, jParams, jMeasurements, jPartition, jObservables);
        }
        
    }
    
}


#endif









