#include "Data.h"
#include "Algebra.h"

namespace data {
    
    
    template <typename Value>
    template <typename Mode>
    Data<Value>::Data(jsx::value const& jParams, Mode) {
        ut::beta = ut::Beta(jParams("beta").real64());         // Bad Bad Bad ...

        filling_ = jsx::at<io::rvec>(jParams("hloc")("filling")); filling_.insert(filling_.begin(), 0);
        if(jParams.is("dyn")) dyn_.reset(new imp::Simple(jParams, jParams("dyn")));
        eig_.reset(new imp::EigenValues<Mode>(jParams, jParams("hloc")("eigen values"), filling(), dyn()));
        ide_.reset(new imp::Operator<Mode, Value>('1', eig()));
        ops_.reset(new imp::Operators<Mode, Value>(jParams, jParams("operators"), eig()));
        hyb_.reset(new bath::Hyb<Value>(jParams, jParams("hybridisation")("matrix"), jParams("hybridisation")("functions")));
    }
    
    template<typename Mode, typename Value>
    void setup_data(jsx::value const& jParams, data::Data<Value>& data)
    {
        cfg::for_each_type<cfg::Worm>::apply(init_worm_data_functor<Mode, Value>(), jParams, data);
    };
    
    template struct Data<double>;
    template struct Data<ut::complex>;
    
    template Data<double>::Data(jsx::value const& jParams, imp::Host);
    template Data<ut::complex>::Data(jsx::value const& jParams, imp::Host);
    
    template void setup_data<imp::Host,double>(jsx::value const& jParams, data::Data<double>& data);
    template void setup_data<imp::Host,ut::complex>(jsx::value const& jParams, data::Data<ut::complex>& data);
    
#ifdef MAKE_GPU_ENABLED
    template Data<double>::Data(jsx::value const& jParams, imp::Device);
    template Data<ut::complex>::Data(jsx::value const& jParams, imp::Device);
#endif
    
    
    
}
