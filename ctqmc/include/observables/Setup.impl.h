#ifndef CTQMC_INCLUDE_OBSERVABLES_SETUP_IMPL_H
#define CTQMC_INCLUDE_OBSERVABLES_SETUP_IMPL_H

#include "Setup.h"

namespace obs {

    
    template<typename Mode, typename Value>
    template<typename W>
    void setup_worm_obs_functor<Mode,Value>::operator()(ut::wrap<W> w, jsx::value const& jParams, data::Data<Value>& data, Observables<Value>& observables) const {
            
            if(jParams.is(W::name())) {
                mpi::cout << "Begin setting up " + W::name() + " observables" << std::endl;
                
                setup_worm_obs<Mode>(jParams, data, observables[get_index<W>::value], w);
                
                mpi::cout << "End setting up " + W::name() + " observables" << std::endl;
            }
            
        }
    
}

#endif
