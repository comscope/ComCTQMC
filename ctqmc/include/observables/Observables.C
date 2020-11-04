#include "Observables.h"

namespace obs {
    
    template<typename Value>
    WormObservables<Value>::WormObservables(std::string worm) : worm_(worm), steps_(0), it_(obs_.end()) {}

    template<typename Value>
    bool WormObservables<Value>::sample(data::Data<Value> const& data, state::State<Value>& state) {
        if(it_ != obs_.end()) return false;
        
        ++steps_;
        
        for(auto const& obs : obs_)
            if(steps_%obs.first == 0)
            {
                sign_ = state.sign();
                it_   = obs_.begin();
                
                return true;
            }
        
        return false;
    }

    template<typename Value>
    bool WormObservables<Value>::cycle(data::Data<Value> const& data, state::State<Value>& state, jsx::value& measurements, imp::itf::Batcher<Value>& batcher) {
        while(it_ != obs_.end())
        {
            if(steps_%it_->first == 0)
                if(!it_->second->sample(sign_, data, state, measurements[worm_], batcher))
                    return false;

            ++it_;
        }
        
        return true;
    }

    template<typename Value>
    void WormObservables<Value>::finalize(data::Data<Value> const& data, jsx::value& measurements) {
        for(auto& obs : obs_)
            obs.second->finalize(data, measurements[worm_]);
    };
    
    template struct WormObservables<double>;
    template struct WormObservables<ut::complex>;

}
