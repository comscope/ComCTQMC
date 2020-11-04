#ifndef CTQMC_INCLUDE_STATE_H
#define CTQMC_INCLUDE_STATE_H

#include <vector>

#include "Utilities.h"
#include "Data.h"
#include "impurity/Product.h"
#include "impurity/Dynamic.h"
#include "impurity/DensityMatrix.h"
#include "impurity/Fact.h"
#include "bath/Bath.h"
#include "config/Expansion.h"
#include "config/Worms.h"


namespace state {
    
    template<typename Mode>
    struct insert_worm_functor {
        
        template<typename Worm, typename Value>
        void operator()(Worm const& worm, data::Data<Value> const& data, state::State<Value>& state);

    };
    
    
    template<typename Value>
    struct State {
        State() = delete;
        template<typename Mode>
        State(jsx::value const& jParams, data::Data<Value> const& data, jsx::value& jConfig, Mode);
        State(State const&) = delete;
        State(State&&) = delete;
        State& operator=(State const&) = delete;
        State& operator=(State&&) = delete;
        ~State() = default;
        
    public:
        
        int& signTimeOrder() {
            return signTimeOrder_;
        };
        cfg::Expansion& expansion() {
            return expansion_;
        };
        cfg::Worm& worm() {
            return worm_;
        };
        imp::itf::Product<Value>& product() {
            return *product_;
        };
        imp::itf::DensityMatrix<Value>& densityMatrix() {
            return *densityMatrix_;
        };
        std::vector<bath::Bath<Value>>& baths() {
            return baths_;
        };
        imp::itf::Dynamic& dyn() {
            return *dyn_;
        };
        imp::Fact<Value>& fact() {
            return fact_;
        };
        
        
        int signTimeOrder() const {
            return signTimeOrder_;
        };
        cfg::Expansion const& expansion() const {
            return expansion_;
        };
        cfg::Worm const& worm() const {
            return worm_;
        };
        imp::itf::Product<Value> const& product() const {
            return *product_;
        };
        imp::itf::DensityMatrix<Value> const& densityMatrix() const {
            return *densityMatrix_;
        };
        std::vector<bath::Bath<Value>> const& baths() const {
            return baths_;
        };
        imp::itf::Dynamic const& dyn() const {
            return *dyn_;
        };
        imp::Fact<Value> const& fact() const {
            return fact_;
        };
        
        
        Value sign() const;
        
        void clean(data::Data<Value> const& data);
        
        jsx::value json() const;
        
    private:
        bool safe_to_load_from_json_;
        int signTimeOrder_;
        cfg::Expansion expansion_;
        cfg::Worm worm_;
        
        std::unique_ptr<imp::itf::Product<Value>> product_;
        std::unique_ptr<imp::itf::DensityMatrix<Value>> densityMatrix_;
        std::vector<bath::Bath<Value>> baths_;
        std::unique_ptr<imp::itf::Dynamic> dyn_;
        imp::Fact<Value> fact_;
    };
    
    
    
    namespace itf {
        
        template<typename Value>
        struct Init {
            virtual bool apply(data::Data<Value> const&, State<Value>&, imp::itf::Batcher<Value>&) = 0;
            virtual ~Init() = default;
        };
        
    }
    
    template<typename Mode, typename Value>
    struct Init : itf::Init<Value> {
        Init() : flag_(ut::Flag::Reject) {};
        Init(Init const&) = delete;
        Init(Init&&) = delete;
        Init& operator=(Init const&) = delete;
        Init& operator=(Init&&) = delete;
        ~Init() = default;
        
        bool apply(data::Data<Value> const& data, State<Value>& state, imp::itf::Batcher<Value>& batcher);
        
    private:
        ut::Flag flag_;
        imp::DensityMatrix<Mode, Value> densityMatrix_;
        
        ut::Flag prepare(data::Data<Value> const& data, State<Value>& state);

        
        ut::Flag decide(State<Value>& state, imp::itf::Batcher<Value>& batcher);
    };

}

#include "State.impl.h"

#endif
