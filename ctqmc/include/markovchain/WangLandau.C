#include "WangLandau.h"

namespace mch {
    
    template<typename Value>
    WangLandau<Value>::WangLandau(jsx::value const& jParams, data::Data<Value> const& data) :
    names_(cfg::Worm::get_names()),
    restart_(jParams.is("restart") and jParams("restart").boolean()),
    flatCriterion_(0.4), lambda_(2.),
    thermalised_(false), totalSteps_(0) {
        steps_.fill(0); eta_.fill(1.);
        
        int index = 0;
        for(auto name : names_) {
            
            if(jParams.is(name)){
                active_.push_back(index);

                if (restart_){
                    if(jParams("measurements").is(name) and jParams("measurements")(name).is("eta")){
                        eta_[index] = jParams("measurements")(name)("eta").real64();
                        steps_[index] = jParams("measurements")(name)("steps").int64() / mpi::number_of_workers();
                        
                    } else {
                        throw std::runtime_error("WangLandau:: restart:: measurement of space " + name + " does not report previous eta.\n No worm space can be added for a `restart' run.");
                    }
                }
            }
            ++index;
        }
        
        if(std::find(active_.begin(), active_.end(), index = get_index<cfg::partition::Worm>::value) == active_.end())
            throw std::runtime_error("mch::WangLandau: partition space not active !");
    };
    
    template<typename Value>
    void WangLandau<Value>::thermalised() {
        if(!thermalised_) {
            if (!restart_){
                normalise_eta();  steps_.fill(0);
            }
            
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
    
    
    template<typename Value>
    void WangLandau<Value>::finalize(jsx::value& measurements) {
        for(auto active : active_){
            measurements(names_[active])["steps"] = steps_[active];
            measurements(names_[active])["eta"] = eta_[active];
        }
    };
    
    template<typename Value>
    jsx::value WangLandau<Value>::etas() {
        jsx::value jEtas;
        
        for(auto active : active_)
            jEtas[names_[active]] = eta_[active];
        
        return jEtas;
    };
    
    template<typename Value>
    void WangLandau<Value>::normalise_eta() {
        double norm = .0;
        for(auto active : active_) norm += eta_[active]*eta_[active];
        
        norm = std::sqrt(norm);
        for(auto active : active_) eta_[active] /= norm;
    }
    
    template struct WangLandau<double>;
    template struct WangLandau<ut::complex>;
    
}

