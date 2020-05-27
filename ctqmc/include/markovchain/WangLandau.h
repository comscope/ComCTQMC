#ifndef CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_H
#define CTQMC_INCLUDE_MARKOVCHAIN_WANGLANDAU_H


#include "../Data.h"
#include "../config/Worms.h"
#include "../../../include/JsonX.h"

namespace mch {
    
    
    template<typename Value>
    struct WangLandau {
        WangLandau() = delete;
        WangLandau(jsx::value const& jParams, data::Data<Value> const& data) :
        names_(cfg::Worm::get_names()),
        flatCriterion_(0.4), lambda_(2.),
        thermalised_(false), totalSteps_(0) {
            steps_.fill(0); eta_.fill(1.);
            //TODO: restart -- get old eta's back and guard against updating them

            int index = 0;
            for(auto name : names_) {
                if(jParams.is(name))  active_.push_back(index);
                ++index;
            }
            
            if(std::find(active_.begin(), active_.end(), index = get_index<cfg::partition::Worm>::value) == active_.end())
                throw std::runtime_error("mch::WangLandau: partition space not active !");
        };
        WangLandau(WangLandau const&) = delete;
        WangLandau(WangLandau&&) = delete;
        WangLandau& operator=(WangLandau const&) = delete;
        WangLandau& operator=(WangLandau&&) = delete;
        ~WangLandau() = default;
        
        void thermalised() {
            if(!thermalised_) {
                normalise_eta();  steps_.fill(0);
            
                //All mp images should have the same eta
                for(auto active : active_) {
                    mpi::reduce<mpi::op::sum>(eta_[active], mpi::master);
                
                    if(mpi::rank() == mpi::master)
                        eta_[active] = eta_[active]/mpi::number_of_workers();
                
                    mpi::bcast(eta_[active], mpi::master);
                }
            
                //Let user know what eta's were chosen
                double const partition = eta_[get_index<cfg::partition::Worm>::value];
                for(auto active : active_) {
                    eta_[active] /= partition;  // Tradition ...
                
                    mpi::cout << names_[active] << " eta = " << eta_[active] << std::endl;
                }

                thermalised_ = true;
            }
        };
        
        template<typename W>
        double eta(W const& w) const {
            return eta_[get_index<W>::value];
        };
        
        template<typename W>
        void inc(W const& w) {
            ++steps_[get_index<W>::value]; ++totalSteps_;

            if(!thermalised_) {
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
        
        void finalize(jsx::value& measurements) {
            for(auto active : active_){
                measurements(names_[active])["steps"] = steps_[active];
                measurements(names_[active])["eta"] = eta_[active];
            }
        };
        
        jsx::value etas() {
            jsx::value jEtas;
            
            for(auto active : active_)
                jEtas[names_[active]] = eta_[active];
            
            return jEtas;
        };
        
    private:
        template<typename W> using get_index = cfg::get_index<W, cfg::Worm>;
        
        std::vector<std::string> const names_;
        double const flatCriterion_;  // Standard deviation at which point lambda is updated
        
        double lambda_;               //Update the volume as  V -> V*lambda
        bool thermalised_;
        std::int64_t totalSteps_;
        
        std::array<double, cfg::Worm::size()> eta_;
        std::array<std::int64_t, cfg::Worm::size()> steps_;
        std::vector<int> active_;
        
        
        void normalise_eta() {
            double norm = .0;
            for(auto active : active_) norm += eta_[active]*eta_[active];
            
            norm = std::sqrt(norm);
            for(auto active : active_) eta_[active] /= norm;
        }
    };
    
}

#endif
