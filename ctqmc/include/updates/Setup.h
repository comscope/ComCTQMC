#ifndef CTQMC_INCLUDE_UPDATES_SETUP_H
#define CTQMC_INCLUDE_UPDATES_SETUP_H

#include <vector>
#include <array>

#include "expansion/Setup.h"

#include "worm/InsertOps.h"
#include "worm/RemoveOps.h"
#include "worm/Reconnect.h"

#include "../markovchain/Update.h"
#include "../markovchain/MarkovChain.h"


namespace upd {
    
    
    template<typename Update, typename Mode, typename Value, typename...Args>
    mch::unique_update_ptr<Value> make_update(double prob, Args&&... args) {
        return mch::unique_update_ptr<Value>(new Generic<Update, Mode, Value>(prob, std::forward<Args>(args)...));
    }
    
    
    template<typename Mode, typename Value>
    void setup_updates(jsx::value const& jParams, data::Data<Value>& data, state::State<Value> const& state, mch::MarkovChain<Value>& markovChain, int const stream);
    
}


#endif
