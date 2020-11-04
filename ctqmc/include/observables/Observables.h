#ifndef CTQMC_INCLUDE_OBSERVABLES_OBSERVABLES_H
#define CTQMC_INCLUDE_OBSERVABLES_OBSERVABLES_H

#include <vector>

#include "Observable.h"
#include "../Data.h"
#include "../State.h"

namespace obs {

    template<typename Value>
    struct WormObservables {
        WormObservables() = delete;
        WormObservables(std::string worm);
        WormObservables(WormObservables const&) = delete;
        WormObservables(WormObservables&&) = default;
        WormObservables& operator=(WormObservables const&) = delete;
        WormObservables& operator=(WormObservables&&) = delete;
        ~WormObservables() = default;
        
        
        template<typename T, typename... Args>
        void add(std::int64_t sweep, std::int64_t store, Args&&... args);
        
        bool sample(data::Data<Value> const& data, state::State<Value>& state);
        
        bool cycle(data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher);
        
        void finalize(data::Data<Value> const& data, jsx::value& measurements);
        
    private:
        using obs_pointer = std::unique_ptr<itf::Observable<Value>>;
        
        std::string  const worm_;
        std::int64_t steps_;
        
        Value sign_;
        std::vector<std::pair<std::int64_t, obs_pointer>> obs_;
        typename std::vector<std::pair<std::int64_t, obs_pointer>>::iterator it_;
    };
    
    
    template<typename Value> using Observables = std::array<std::unique_ptr<WormObservables<Value>>, cfg::Worm::size()>;

}

#include "Observables.impl.h"

#endif
