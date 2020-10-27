#include "State.h"
#include "Algebra.h"

namespace state {
    
    template<typename Value>
    template<typename Mode>
    State<Value>::State(jsx::value const& jParams, data::Data<Value> const& data, jsx::value& jConfig, Mode) :
    safe_to_load_from_json_(!jConfig.is("worm") or jParams.is(jConfig("worm")("name").string()) ? true : false),
    signTimeOrder_(1),
    expansion_( safe_to_load_from_json_ and jConfig.is("expansion") ? jConfig("expansion") : jsx::null_t(), data.ops().flavors()),
    worm_( safe_to_load_from_json_ and jConfig.is("worm") ? jConfig("worm") : jsx::object_t{{"name", "partition"}, {"entry", jsx::null_t()}}),
    product_(new imp::Product<Mode, Value>(jParams, data.eig(), data.ide(), data.ops())),
    densityMatrix_(new imp::DensityMatrix<Mode, Value>()),
    baths_(data.hyb().blocks().size()),
    dyn_(data.dyn() != nullptr ? new imp::Dynamic(*data.dyn()) : new imp::itf::Dynamic()) {
        int index = 0;
        for(auto const& block : data.hyb().blocks()) {
            for(auto const& flavor : block.flavorsL())
                for(auto const& key : expansion_.at(flavor))
                    baths_[index].insertL(key, flavor);
            
            for(auto const& flavor : block.flavorsR())
                for(auto const& key : expansion_.at(flavor))
                    baths_[index].insertR(key, flavor);
            
            baths_[index++].clean(data.hyb());
        }
        
        for(auto const& bath : baths()) {
            auto const& opsL = bath.opsL();
            auto const& opsR = bath.opsR();
            
            for(std::size_t i = 0; i < opsL.size(); ++i) {
                if(!product().insert(opsR[i].key(), opsR[i].flavor()))
                    throw std::runtime_error("state::State::constructor: key in config appears twice.");
                if(!product().insert(opsL[i].key(), opsL[i].flavor()))
                    throw std::runtime_error("state::State::constructor: key in config appears twice.");
                
                dyn().insert(opsR[i].key(), opsR[i].flavor());
                dyn().insert(opsL[i].key(), opsL[i].flavor());
            }
        }

        cfg::apply(worm(), insert_worm_functor<Mode>(), data, *this);


        dyn().ratio();
        dyn().accept();
        
        fact().accept();
    };
    
    
    
    template<typename Value>
    Value State<Value>::sign() const {
        Value sign = static_cast<double>(signTimeOrder())*densityMatrix().sign()*fact().sign();  //!!!!!!!!!!!!!!!!
        for(auto const& bath : baths()) sign *= bath.sign();
        return sign;
    };
    
    template<typename Value>
    void State<Value>::clean(data::Data<Value> const& data) {
        dyn().clean();
        for(auto& bath : baths())
            bath.clean(data.hyb());
    };
    
    template<typename Value>
    jsx::value State<Value>::json() const {
        return jsx::object_t{
            {"expansion", expansion().json()},
            {"worm", worm().json()}
        };
    };
    
    template struct State<double>;
    template struct State<ut::complex>;
    
    template State<double>::State(jsx::value const& jParams, data::Data<double> const& data, jsx::value& jConfig, imp::Host);
    template State<ut::complex>::State(jsx::value const& jParams, data::Data<ut::complex> const& data, jsx::value& jConfig, imp::Host);
    
#ifdef MAKE_GPU_ENABLED
    template State<double>::State(jsx::value const& jParams, data::Data<double> const& data, jsx::value& jConfig, imp::Device);
    template State<ut::complex>::State(jsx::value const& jParams, data::Data<ut::complex> const& data, jsx::value& jConfig, imp::Device);
#endif
    
    template<typename Mode, typename Value>
    bool Init<Mode,Value>::apply(data::Data<Value> const& data, State<Value>& state, imp::itf::Batcher<Value>& batcher) {
        if(flag_ != ut::Flag::Pending) flag_ = prepare(data, state);
        if(flag_ == ut::Flag::Pending) flag_ = decide(state, batcher);
        return flag_ != ut::Flag::Pending;
    };
    
    template<typename Mode, typename Value>
    ut::Flag Init<Mode,Value>::prepare(data::Data<Value> const& data, State<Value>& state) {
        densityMatrix_ = imp::DensityMatrix<Mode, Value>(state.product(), data.eig());
        if(ut::Flag::Pending != densityMatrix_.surviving(data.eig()))
            throw std::runtime_error("state::Init: initial trace is zero");
        return ut::Flag::Pending;
    };
    
    template<typename Mode, typename Value>
    ut::Flag Init<Mode,Value>::decide(State<Value>& state, imp::itf::Batcher<Value>& batcher) {
        if(ut::Flag::Pending == densityMatrix_.decide(.0, state.product(), batcher)) return ut::Flag::Pending;
        imp::get<Mode>(state.densityMatrix()) = std::move(densityMatrix_);
        state.signTimeOrder() *= state.product().accept();
        
        return ut::Flag::Accept;
    };
    
    template struct Init<imp::Host,double>;
    template struct Init<imp::Host,ut::complex>;
    
    
#ifdef MAKE_GPU_ENABLED
    template struct Init<imp::Device,double>;
    template struct Init<imp::Device,ut::complex>;
#endif
    
}
