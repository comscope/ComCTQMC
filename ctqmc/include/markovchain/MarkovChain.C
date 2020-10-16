#include "../Algebra.h"
#include "MarkovChain.h"

namespace mch {
    
    std::int64_t select_seed(jsx::value const& jParams, std::int64_t const ID){
        
        auto inc = (jParams.is("seed increment") ? jParams("seed increment").int64() : 857)*ID;
        std::int64_t seed = 0;
        
        if (jParams.is("restart") and jParams("restart").boolean()){
            seed = static_cast<int64_t>(time(0));
            
        } else {
            seed = jParams.is("seed") ? jParams("seed").int64() : 41085;
            
        }
        
        return seed + inc;
    }
    
    
    template <typename Value>
    template<typename Mode>
    MarkovChain<Value>::MarkovChain(jsx::value const& jParams, std::int64_t ID, Mode) :
        clean_(jParams.is("clean") ? jParams("clean").int64() : 10000), steps_(0),
        init_(new state::Init<Mode, Value>()),
        urng_(ut::Engine(select_seed(jParams,ID)), ut::UniformDistribution(.0, 1.)),
        update_(nullptr) {
        };
    
    
    template <typename Value>
    void MarkovChain<Value>::add(std::unique_ptr<itf::Update<Value>> update) {
        if(update->origin() != update->target())
            throw std::runtime_error("mc::MarkovChain::add");
        
        add_update(std::move(update));
    };
    
    template <typename Value>
    void MarkovChain<Value>::add(std::unique_ptr<itf::Update<Value>> updateAB, std::unique_ptr<itf::Update<Value>> updateBA) {
        if(updateAB->origin() != updateBA->target() || updateBA->origin() != updateAB->target())
            throw std::runtime_error("mc::MarkovChain::add");
        
        updateAB->ratioChoose_ = updateBA->probChoose()/updateAB->probChoose();  //not yet properly normalized
        updateBA->ratioChoose_ = updateAB->probChoose()/updateBA->probChoose();  //not yet properly normalized
        
        add_update(std::move(updateAB));
        add_update(std::move(updateBA));
    };
    
    template <typename Value>
    void MarkovChain<Value>::finalize(state::State<Value> const& state) {
        for(auto& updates : allUpdates_)
            for(auto& update : updates)
                update->ratioChoose_ *= allDistrs_[update->origin()].back()/allDistrs_[update->target()].back();  // now they are properly normalized
        
        choose_update(state);
    };
    
    template <typename Value>
    bool MarkovChain<Value>::init(data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
        return init_->apply(data, state, batcher);
    };
    
    template <typename Value>
    bool MarkovChain<Value>::cycle(mch::WangLandau<Value>& wangLandau,  data::Data<Value> const& data, state::State<Value>& state, imp::itf::Batcher<Value>& batcher) {
        if(!update_->apply(urn_, wangLandau, data, state, urng_, batcher)) return false;
        
        choose_update(state); if(++steps_ % clean_ == 0) state.clean(data);
        
        return true;
    };
    
    template <typename Value>
    void MarkovChain<Value>::add_update(std::unique_ptr<itf::Update<Value>> update) {
        auto& distr = allDistrs_[update->origin()];
        auto& updates = allUpdates_[update->origin()];
        
        distr.push_back((distr.size() ? distr.back() : .0) + update->probChoose());
        updates.push_back(std::move(update));
    };
    
    template <typename Value>
    void MarkovChain<Value>::choose_update(state::State<Value> const& state) {
        auto& distr = allDistrs_[state.worm().index()];
        auto& updates = allUpdates_[state.worm().index()];
        
        double const urn = urng_()*distr.back();
        auto const it = std::upper_bound(distr.begin(), distr.end(), urn);
        double const low = it != distr.begin() ? *(it - 1) : .0;
        
        urn_ = (urn - low)/(*it - low);
        update_ = updates[it - distr.begin()].get();
    };
        
    template struct MarkovChain<double>;
    template struct MarkovChain<ut::complex>;
    
    template MarkovChain<double>::MarkovChain(jsx::value const& jParams, std::int64_t ID, imp::Host);
    template MarkovChain<ut::complex>::MarkovChain(jsx::value const& jParams, std::int64_t ID, imp::Host);
    
#ifdef MAKE_GPU_ENABLED
    template MarkovChain<double>::MarkovChain(jsx::value const& jParams, std::int64_t ID, imp::Device);
    template MarkovChain<ut::complex>::MarkovChain(jsx::value const& jParams, std::int64_t ID, imp::Device);
#endif
    
}

