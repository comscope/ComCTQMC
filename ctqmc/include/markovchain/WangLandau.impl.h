#ifndef CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_IMPL_H
#define CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_IMPL_H


#include "WangLandau.h"

namespace mch {

template<typename Value>
template<typename W>
void WangLandau<Value>::inc(W const& w) {
    if (!restart_ or thermalised_){
        ++steps_[get_index<W>::value]; ++totalSteps_;
    }
    
    if(!thermalised_ and !restart_) {
        //Update eta, avoid underflows ...
        if(eta_[get_index<W>::value] > 1.e-25)
            eta_[get_index<W>::value] /= lambda_;
        else {
            normalise_eta();
            
            if(eta_[get_index<W>::value] > 1.e-25)
                eta_[get_index<W>::value] /= lambda_;
            else {
                lambda_ = std::sqrt(lambda_);  eta_.fill(1.);  steps_.fill(0);  totalSteps_ = 0;  return; // ... and restart with smaller lambda if eta's get crazy
            }
        }
        
        //Look at the current histogram ...
        std::int64_t min = std::numeric_limits<std::int64_t>::max(), max = 0;
        for(auto active : active_) {
            min = std::min(steps_[active], min);  max = std::max(steps_[active], max);
        }
        
        // ... is it sufficiently flat?  if so, update lambda and reset histogram
        if(max - min < totalSteps_*flatCriterion_/active_.size()) {
            lambda_ = std::sqrt(lambda_);  steps_.fill(0);  totalSteps_ = 0;
        }
    }
    
};
    
}

#endif
