#ifndef CTQMC_INCLUDE_DATA_IMPL_H
#define CTQMC_INCLUDE_DATA_IMPL_H

#include "Data.h"

namespace data {
    
    template<typename Mode, typename Value>
    template<typename W>
    void init_worm_data_functor<Mode, Value>::operator()(ut::wrap<W> w, jsx::value const& jParams, data::Data<Value>& data) const {
        if(jParams.is(W::name()))
            cfg::init_worm_data<Mode>(w, jParams, data);
    }


}

#endif
