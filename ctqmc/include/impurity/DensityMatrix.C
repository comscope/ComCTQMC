#include "../Algebra.h"
#include "DensityMatrix.h"

namespace imp {
    
    //-----------------------------------------------------------------DENSITYMATRIX---------------------------------------------------------------------------------
    
    
    template<typename Mode, typename Value>
    DensityMatrix<Mode, Value>::DensityMatrix(itf::Product<Value>& product, itf::EigenValues const& eig) :
        level_(product.height()),
        Z_(.0), z_(eig.sectorNumber() + 1) {
            std::vector<int> source(eig.sectorNumber()); std::iota(source.begin(), source.end(), 1);
            std::vector<int> target = get<Mode>(product).map(get<Mode>(product).first(), level_, source);

            for(int i = 0; i < eig.sectorNumber(); ++i)
                if(target[i] == source[i])
                    bounds_.push_back({source[i], get<Mode>(product).first()->op(level_)->map(source[i]).norm, ut::Zahl<double>()});
        };
    
    template<typename Mode, typename Value>
    ut::Flag DensityMatrix<Mode, Value>::surviving(itf::EigenValues const& eig) {
        if(bounds_.size() == 0) return ut::Flag::Reject;

        for(auto it = bounds_.begin(); it != bounds_.end(); ++it) it->ln += get<Mode>(eig).at(it->sec).ln_dim();

        std::sort(bounds_.begin(), bounds_.end(), &compare);
        
        for(auto it = bounds_.begin(); it != bounds_.end(); ++it) it->value = ut::exp(it->ln);
        for(auto it = bounds_.end() - 1; it != bounds_.begin(); --it) (it - 1)->value += it->value;
        bounds_.push_back({0, .0, ut::Zahl<double>(.0)});

        bound_ = bounds_.begin(); return ut::Flag::Pending;
    };
    
    template<typename Mode, typename Value>
    ut::Flag DensityMatrix<Mode, Value>::decide(ut::Zahl<double> const& thresh, itf::Product<Value>& product, itf::Batcher<Value>& batcher) {
        if(ut::abs(Z_) + bound_->value <= ut::abs(thresh)) {
            return ut::Flag::Reject;
        } else if(bound_->value <= std::numeric_limits<double>::epsilon()*ut::abs(Z_)) {
            op_ = get<Mode>(product).first()->get_op(level_);
            return ut::Flag::Accept;
        }

        sectors_.push_back(bound_->sec);
        get<Mode>(product).multiply(get<Mode>(product).first(), level_, bound_->sec, batcher);
        trace(&z_[bound_->sec], &Z_, get<Mode>(product).first()->op(level_)->mat(bound_->sec), batcher);
        
        ++bound_; return ut::Flag::Pending;
    };
    
    template<typename Mode, typename Value>
    Value DensityMatrix<Mode, Value>::sign() const { return Z_.mantissa()/std::abs(Z_.mantissa());};
    
    template struct DensityMatrix<imp::Host, double>;
    template struct DensityMatrix<imp::Host, ut::complex>;
    
#ifdef MAKE_GPU_ENABLED
    template struct DensityMatrix<imp::Device, double>;
    template struct DensityMatrix<imp::Device, ut::complex>;
#endif
    
}
