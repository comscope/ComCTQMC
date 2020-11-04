#ifndef CTQMC_INCLUDE_OBSERVABLES_OBSERVABLES_IMPL_H
#define CTQMC_INCLUDE_OBSERVABLES_OBSERVABLES_IMPL_H

#include "Observables.h"

namespace obs {

    template<typename Value>
    template<typename T, typename... Args>
    void WormObservables<Value>::add(std::int64_t sweep, std::int64_t store, Args&&... args) {
            obs_.emplace_back(sweep, obs_pointer(new T(store, std::forward<Args>(args)...)));
            
            it_ = obs_.end();
    };
        
}

#endif
