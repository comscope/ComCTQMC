#ifndef INCLUDE_PARAMETERS_INITIALIZE_IMPL_H
#define INCLUDE_PARAMETERS_INITIALIZE_IMPL_H

#include <complex>

#include "Defaults.h"
#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../../ctqmc/include/config/Worms.h"

namespace params {

        template<typename W>
        void worm_defaults_functor::operator()(ut::wrap<W> w, jsx::value & jParams, AllDefaults const& defaults) const {
            
            if( W::name() != cfg::partition::Worm::name() && jParams.is(W::name()) ) {
                
                if (cfg::worm_size<W>::value == 2){
                    for (auto const& entry : defaults.oneTimeWormDefaults.get().object()){
                        if (!jParams(W::name()).is(entry.first))
                            jParams[W::name()][entry.first] = entry.second;
                        else check_type_against_default(W::name() + " -> " + entry.first, jParams[W::name()][entry.first], entry.second);
                    }
                } else {
                    for (auto const& entry : defaults.multiTimeWormDefaults.get().object()){
                        if (!jParams(W::name()).is(entry.first))
                            jParams[W::name()][entry.first] = entry.second;
                        else check_type_against_default(W::name() + " -> " + entry.first, jParams[W::name()][entry.first], entry.second);
                    }
                }
            }
        }


        template<typename W>
           void worm_cutoffs_functor::operator()(ut::wrap<W> w, jsx::value & jParams, int const fermion_cutoff, int const boson_cutoff) const {
                
               if( W::name() != cfg::partition::Worm::name() && jParams.is(W::name()) ) {
                   
                   auto & jWorm = jParams[W::name()];
                   
                   if (jWorm.is("fermion cutoff") and jWorm("basis").string() == "matsubara" ) jWorm["fermion cutoff"] = fermion_cutoff;
                   if (jWorm.is("matsubara cutoff")) jWorm["matsubara cutoff"] = fermion_cutoff;
                   
                   if (jWorm.is("cutoff")) jWorm["cutoff"] = boson_cutoff;
                   if (jWorm.is("boson cutoff")) jWorm["boson cutoff"] = boson_cutoff;
                       
                
               }
               
           }

}

#endif
