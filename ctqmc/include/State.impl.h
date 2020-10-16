#ifndef CTQMC_INCLUDE_STATE_IMPL_H
#define CTQMC_INCLUDE_STATE_IMPL_H

#include "State.h"

namespace state {

    template<typename Mode>
    template<typename Worm, typename Value>
    void insert_worm_functor<Mode>::operator()(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state) {
        if(!cfg::insert_worm<Mode>(worm, data, state))
            throw std::runtime_error("state::insert_worm_functor: key appears twice in " + Worm::name() + " worm");
    };
    
}

#endif
