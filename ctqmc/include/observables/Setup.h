#ifndef CTQMC_INCLUDE_OBSERVABLES_SETUP_H
#define CTQMC_INCLUDE_OBSERVABLES_SETUP_H

#include <vector>

#include "partition/Setup.h"
#include "worm/Setup.h"


namespace obs {

    
    template<typename Mode, typename Value>
    struct setup_worm_obs_functor {
        template<typename W> using get_index = cfg::get_index<W, cfg::Worm>;
        
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value const& jParams, data::Data<Value>& data, Observables<Value>& observables) const;
    };

    template<typename Mode, typename Value>
    void setup_obs(jsx::value const& jParams, data::Data<Value>& data, Observables<Value>& observables);
}

#include "Setup.impl.h"

#endif
