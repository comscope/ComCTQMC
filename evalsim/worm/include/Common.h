#ifndef EVALSIM_WORM_INCLUDE_COMMON
#define EVALSIM_WORM_INCLUDE_COMMON


#include "../../../include/JsonX.h"
#include "../../../ctqmc/include/config/Worms.h"
#include "functions/Functions.h"
#include "functions/Measurements.h"
#include "functions/Utilities.h"



namespace evalsim {
    
    namespace worm {

        namespace func{
            
            template <typename Value>
            struct BosonFrequencies {
                
                BosonFrequencies() = delete;
                BosonFrequencies(jsx::value const& jWorm) :
                pos_and_neg_freq_(std::is_same<Value,double>::value ? 1 : 2),
                nMatGB_(pos_and_neg_freq_*jWorm("cutoff").int64()-pos_and_neg_freq_+1),
                omega_b_(nMatGB_,true,!std::is_same<Value,double>::value){}
                
                int nMatGB() const {return nMatGB_;}
                
                OmegaMap omega_b() const {return omega_b_;}
                
            private:
                int const pos_and_neg_freq_, nMatGB_;
                OmegaMap const omega_b_;
            };
        
            template <typename Value>
            struct Frequencies {
                
                Frequencies() = delete;
                Frequencies(jsx::value const& jWorm) :
                pos_and_neg_freq_(std::is_same<Value,double>::value ? 1 : 2),
                nMatGB_(pos_and_neg_freq_*jWorm("boson cutoff").int64()-pos_and_neg_freq_+1),
                nMatGF_(2*(jWorm("basis").string() == "matsubara" ? jWorm("fermion cutoff").int64() : ( jWorm.is("matsubara cutoff") ? jWorm("matsubara cutoff").int64() : 50 ))),
                omega_b_(nMatGB_,true,!std::is_same<Value,double>::value),
                omega_f_(nMatGF_,false,true){}
                
                int nMatGB() const {return nMatGB_;}
                int nMatGF() const {return nMatGF_;}
                
                OmegaMap omega_b() const {return omega_b_;}
                OmegaMap omega_f() const {return omega_f_;}
                
            private:
                int const pos_and_neg_freq_, nMatGB_, nMatGF_;
                OmegaMap const omega_b_, omega_f_;
            };
        
        }
    }
}

#endif



