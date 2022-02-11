#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_GLOBALSPINFLIP_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_GLOBALSPINFLIP_H

namespace upd {
    
namespace expansion {
        

    template<typename Worm>
    struct SpinFlip {
        SpinFlip() = delete;
        template<typename Value>
        SpinFlip(jsx::value const& jParams, data::Data<Value> const& data) : flipMap_(data.hyb().flavors()), prob_(jParams("spin flip").real64()) {
            
            int const N = data.hyb().flavors();
            int const n = N/2;
            
            int bath = 0;
            for(auto const& block : data.hyb().blocks()){
                
                for(auto const& flavorL : block.flavorsL()){
                    
                    flipMap_[flavorL].bath = bath;
                    flipMap_[flavorL].flavor = flavorL;
                    flipMap_[flavorL].flipped_flavor = (flavorL + n) % N; //Doesn't work for coupled basis
                    
                    int flipped_bath = 0;
                    for(auto const& flipped_block : data.hyb().blocks()){
                        for(auto const& flipped_flavor : flipped_block.flavorsL()){
                            if (flipped_flavor == flipMap_[flavorL].flipped_flavor){
                                flipMap_[flavorL].flipped_bath = flipped_bath;
                            }
                        }
                    flipped_bath++;
                    }
                }
                    
                for(auto const& flavorR : block.flavorsR()){
                    
                    flipMap_[flavorR].bath = bath;
                    flipMap_[flavorR].flavor = flavorR;
                    flipMap_[flavorR].flipped_flavor = (flavorR + n) % N; //Doesn't work for coupled basis
                    
                    int flipped_bath = 0;
                    for(auto const& flipped_block : data.hyb().blocks()){
                        for(auto const& flipped_flavor : flipped_block.flavorsR()){
                            if (flipped_flavor == flipMap_[flavorR].flipped_flavor){
                                flipMap_[flavorR].flipped_bath = flipped_bath;
                            }
                        }
                        flipped_bath++;
                    }
                }
                bath++;
            }
            
        }
        
        jsx::value json() const {
            return jsx::null_t();
        };
        
        
        template<typename Mode, typename Value>
        bool impurity(data::Data<Value> const& data, state::State<Value>& state) const {
            
            for (auto const& flip : flipMap_){
                auto & ops = state.expansion()[flip.flavor];
                for (auto & op : ops){
                    
                    state.product().erase(op);
                    if(!state.product().insert(op, flip.flipped_flavor)) return false;
                    
                    state.dyn().erase(op, flip.flavor);
                    state.dyn().insert(op, flip.flipped_flavor);
                    
                }
            }
            
            return true;
        };
        
        template<typename Value>
        void bath(data::Data<Value> const& data, std::vector<bath::Bath<Value>>& baths) const {
            
            for (int i_flip = 0; i_flip < flipMap_.size(); i_flip+=2){
                auto const& opsL = baths[i_flip].opsL();
                auto const& opsR = baths[i_flip+1].opsR();
                
                for (int i_op = 0; i_op < opsL.size(); i_op++){
                    
                    baths[flipMap_[i_flip].bath].add(data.hyb(), bath::Erase{opsL[i_op].key(), opsR[i_op].key()});
                    baths[flipMap_[i_flip].flipped_bath].add(data.hyb(), bath::Insert{
                        opsL[i_op].key(), flipMap_[i_flip].flipped_flavor,
                        opsR[i_op].key(), flipMap_[i_flip+1].flipped_flavor} );
                }
            }
            
        };
        
        template<typename Value>
        void accept(data::Data<Value> const& data, state::State<Value>& state) const {
                        
            for (auto const& flip : flipMap_){
                auto & ops = state.expansion()[flip.flavor];
                for (auto & op : ops){
                    state.expansion()[flip.flavor].erase(op);
                    state.expansion()[flip.flipped_flavor].insert(op);
                }
            }
        };
        
        template<typename Value>
        void reject(data::Data<Value> const& data, state::State<Value>& state) {
        };
        
        template<typename Value>
        bool propose(double const urn, data::Data<Value> const& data, state::State<Value>& state, ut::UniformRng& urng) { return true; };
        
        Worm const& target(Worm const& worm) const {
            return worm;
        };
        
        template<typename Value>
        double ratio(data::Data<Value> const& data, state::State<Value>& state) const {
            return 1.;
        };
        
        protected:
            struct FlipMapElement {
                FlipMapElement() = default;
            
                int flavor, bath, flipped_flavor, flipped_bath;
            };
                            
            std::vector<ut::KeyType> keys;
            double prob_;
            std::vector<FlipMapElement> flipMap_;
        
    };

}
        
}



#endif
