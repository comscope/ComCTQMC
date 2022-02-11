#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_OpShift_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_OpShift_H

#include "../include/Generic.h"
#include "../../bath/algebra/Exchange.h"

namespace upd {

namespace expansion {

    struct ShiftBase {
        ShiftBase() = delete;
        template<typename Value>
        ShiftBase(jsx::value const& jParams, data::Data<Value> const& data) {
            int bath = 0;
            for(auto const& block : data.hyb().blocks()) {
                for(auto const& flavorL : block.flavorsL())
                    for(auto const& flavorR : block.flavorsR()) {
                        prob_.push_back((prob_.size() ? prob_.back() : .0) + 1.);
                        defs_.push_back({flavorL, flavorR, bath});
                    }
                ++bath;
            }
        }
        ~ShiftBase() = default;
        
        jsx::value json() const {
            return jsx::null_t();
        };
        
        template<typename Value>
        double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
            return 1.;
        };
        
        template<typename Value>
        void reject(data::Data<Value> const& data, state::State<Value>& state) {
        };
        
        template<typename Value>
        bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) {
            def_ = defs_[std::upper_bound(prob_.begin(), prob_.end(), urn*prob_.back()) - prob_.begin()];
            
            auto const& opsL = state.expansion()[def_.flavorL];
            auto const& opsR = state.expansion()[def_.flavorR];
            
            if(!(opsL.size() && opsR.size())) return false;
            
            ut::KeyType keyHigh;
            if(opsL.size() + opsR.size() > 2) {
                std::size_t index = urng()*(opsL.size() + opsR.size());
                ut::KeyType const keyLow = index < opsL.size() ? opsL[index] : opsR[index - opsL.size()];
                
                auto itL = std::upper_bound(opsL.begin(), opsL.end(), keyLow); ut::KeyType shiftL = 0;
                if(itL != opsL.end()) def_.keyL = *itL; else { def_.keyL = *(itL = opsL.begin()); shiftL = ut::KeyMax;};
                
                auto itR = std::upper_bound(opsR.begin(), opsR.end(), keyLow); ut::KeyType shiftR = 0;
                if(itR != opsR.end()) def_.keyR = *itR; else { def_.keyR = *(itR = opsR.begin()); shiftR = ut::KeyMax;};
                
                keyHigh = ++itL != opsL.end() ? *itL + shiftL : *opsL.begin() + ut::KeyMax;
                keyHigh = std::min(++itR != opsR.end() ? *itR + shiftR : *opsR.begin() + ut::KeyMax, keyHigh);
                
                if(std::max(def_.keyL + shiftL, def_.keyR + shiftR) < keyHigh) keyDiff_ = keyHigh - keyLow; else return 0;
                
            } else {
                def_.keyL = *opsL.begin();
                def_.keyR = *opsR.begin();
                
                def_.newKey = ut::KeyMax*urng();
                
                keyHigh = keyDiff_ = ut::KeyMax;
            }
            
            def_.newKey = ut::cyclic(keyHigh - urng()*keyDiff_);
            if (def_.newKey > ut::KeyMax) def_.newKey -= ut::KeyMax;
            
            return true;
        };
        
    protected:
        struct Def {
            Def() = default;
            Def(int flavorL, int flavorR, int bath) : flavorL(flavorL), flavorR(flavorR), bath(bath) {};
            int flavorL, flavorR, bath;
            ut::KeyType keyL, keyR, newKey;
        };
        
        std::vector<double> prob_;
        std::vector<Def> defs_;
        Def def_;
        ut::KeyType keyDiff_;
    };

    template <typename Worm>
    struct RShift : ShiftBase {
        template<typename Value>
        RShift(jsx::value const& jParams, data::Data<Value> const& data) : ShiftBase(jParams, data) {};
        
        template<typename Mode, typename Value>
        bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
            
            state.product().erase(def_.keyR);
            if(!state.product().insert(def_.newKey, def_.flavorR)) return false;
            
            state.dyn().erase(def_.keyR, def_.flavorR);
            state.dyn().insert(def_.newKey, def_.flavorR);
            
            return true;
        };
        
        template<typename Value>
        void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
            
            baths[def_.bath].add(data.hyb(), bath::Exchange{def_.keyR, def_.newKey, def_.flavorR});
            
        };
        
        template<typename Value>
        void accept(data::Data<Value> const& data, state::State<Value>& state) const {
            state.expansion()[def_.flavorR].erase(def_.keyR);
            state.expansion()[def_.flavorR].insert(def_.newKey);
        };
               
        Worm const& target(Worm const& worm) const {
            return worm;
        };
    };

    template <typename Worm>
    struct LShift : ShiftBase {
        template<typename Value>
        LShift(jsx::value const& jParams, data::Data<Value> const& data) : ShiftBase(jParams, data) {};
        
        template<typename Mode, typename Value>
        bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
            
            state.product().erase(def_.keyL);
            if(!state.product().insert(def_.newKey, def_.flavorL)) return false;
            
            state.dyn().erase(def_.keyL, def_.flavorL);
            state.dyn().insert(def_.newKey, def_.flavorL);
            
            return true;
        };
        
        template<typename Value>
        void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
            
            baths[def_.bath].add(data.hyb(), bath::Exchange{def_.keyL, def_.newKey, def_.flavorL});
            
        };
        
        template<typename Value>
        void accept(data::Data<Value> const& data, state::State<Value>& state) const {
            state.expansion()[def_.flavorL].erase(def_.keyL);
            state.expansion()[def_.flavorL].insert(def_.newKey);
        };
        
        Worm const& target(Worm const& worm) const {
            return worm;
        };
    };

}

}

#endif
