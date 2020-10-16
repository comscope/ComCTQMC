#ifndef CTQMC_INCLUDE_MARKOVCHAIN_MARKOVCHAIN_H
#define CTQMC_INCLUDE_MARKOVCHAIN_MARKOVCHAIN_H

#include <vector>
#include <array>
#include <ctime>

#include "Update.h"
#include "../Utilities.h"
#include "../Data.h"
#include "../State.h"


namespace mch {
    
    std::int64_t select_seed(jsx::value const& jParams, std::int64_t const ID);

    template<typename Value>
    struct MarkovChain {
        MarkovChain() = delete;
        template<typename Mode>
        MarkovChain(jsx::value const& jParams, std::int64_t ID, Mode);
        MarkovChain(MarkovChain const&) = delete;
        MarkovChain(MarkovChain&&) = delete;
        MarkovChain& operator=(MarkovChain const&) = delete;
        MarkovChain& operator=(MarkovChain&&) = delete;
        ~MarkovChain() = default;
        
        
        void add(std::unique_ptr<itf::Update<Value>> update);
        
        void add(std::unique_ptr<itf::Update<Value>> updateAB, std::unique_ptr<itf::Update<Value>> updateBA);
        
        void finalize(state::State<Value> const& state);
        
        bool init(data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher);
        
        bool cycle(mch::WangLandau<Value>& wangLandau,  data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher);
        
    private:
        std::int64_t const clean_;
        std::int64_t steps_;
        double urn_;
        
        std::unique_ptr<state::itf::Init<Value>> init_;
        
        ut::UniformRng urng_;
        std::array<std::vector<double>, cfg::Worm::size()> allDistrs_;
        std::array<std::vector<std::unique_ptr<itf::Update<Value>>>, cfg::Worm::size()> allUpdates_;
        itf::Update<Value>* update_;

        
        void add_update(std::unique_ptr<itf::Update<Value>> update);
        
        void choose_update(state::State<Value> const& state);
    };
    
}

#include "MarkovChain.impl.h"

#endif
