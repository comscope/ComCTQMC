#ifndef CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_H
#define CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_H


#include "../Data.h"
#include "../config/Worms.h"
#include "../../../include/JsonX.h"

namespace mch {
    
    
    template<typename Value>
    struct WangLandau {
        WangLandau() = delete;
        WangLandau(jsx::value const& jParams, data::Data<Value> const& data);
        WangLandau(WangLandau const&) = delete;
        WangLandau(WangLandau&&) = delete;
        WangLandau& operator=(WangLandau const&) = delete;
        WangLandau& operator=(WangLandau&&) = delete;
        ~WangLandau() = default;
        
        void thermalised();
        
        template<typename W>
        double eta(W const& w) const {
            return eta_[get_index<W>::value];
        };
        
        template<typename W>
        void inc(W const& w);
        
        void finalize(jsx::value& measurements);
        
        jsx::value etas();
        
    private:
        template<typename W> using get_index = cfg::get_index<W, cfg::Worm>;
        
        std::vector<std::string> const names_;
        bool const restart_; //is the job a "restart" job -- in which case we don't recompute eta's
        double const flatCriterion_;  // Standard deviation at which point lambda is updated
        
        double lambda_;               //Update the volume as  V -> V*lambda
        bool thermalised_;
        std::int64_t totalSteps_;
        
        std::array<double, cfg::Worm::size()> eta_;
        std::array<std::int64_t, cfg::Worm::size()> steps_;
        std::vector<int> active_;
        
        
        void normalise_eta();
    };
    
}

#include "WangLandau.impl.h"

#endif
