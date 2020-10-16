
#include "../Algebra.h"
#include "Setup.h"


namespace obs {

    template<typename Mode, typename Value>
    void setup_obs(jsx::value const& jParams, data::Data<Value>& data, Observables<Value>& observables) {

        cfg::for_each_type<cfg::Worm>::apply(setup_worm_obs_functor<Mode, Value>(), jParams, data, observables);
        
    }
    
    template void setup_obs<imp::Host, double>(jsx::value const& jParams, data::Data<double>& data, Observables<double>& observables);
    template void setup_obs<imp::Host, ut::complex>(jsx::value const& jParams, data::Data<ut::complex>& data, Observables<ut::complex>& observables);
    
#ifdef MAKE_GPU_ENABLED
    
    template void setup_obs<imp::Device, double>(jsx::value const& jParams, data::Data<double>& data, Observables<double>& observables);
    template void setup_obs<imp::Device, ut::complex>(jsx::value const& jParams, data::Data<ut::complex>& data, Observables<ut::complex>& observables);
    
#endif

}
