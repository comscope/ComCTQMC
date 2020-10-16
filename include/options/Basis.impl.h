#ifndef INCLUDE_OPTIONS_BASIS_IMPL_H
#define INCLUDE_OPTIONS_BASIS_IMPL_H

#include "Basis.h"

namespace opt {

template <typename Value>
template <typename Type>
Basis<Value>::Basis(Type const& type, jsx::value const& jTransformation) :
n_(type.n()), transformation_(2*n_, jTransformation) {
    
    u_.resize(transformation_.I()*transformation_.J(), .0);
    for(int f_new = 0; f_new < transformation_.I(); ++f_new)
    for(int f_old = 0; f_old < transformation_.J(); ++f_old)
    for(int m = 0; m < n_; ++m)
    for(int s = 0; s < 2; ++s)
    u_[f_new*2*n_ + 2*m + s] += transformation_(f_new, f_old)*type(f_old, m, s);
    
    for(auto const& qn : type.qns()) {
        auto const& qn_old = qn.second; io::rvec qn_new;
        
        for(int f_new = 0; f_new < transformation_.I(); ++f_new) {
            std::set<double> count;
            
            for(int f_old = 0; f_old < transformation_.J(); ++f_old)
            if(transformation_(f_new, f_old) != .0) count.insert(qn_old.at(f_old));
            
            if(count.size() == 1) qn_new.push_back(*count.begin());
        }
        
        if(qn_new.size() == static_cast<std::size_t>(transformation_.I())) qns_[qn.first] = qn_new;
    }
}

};

#endif //OPTIONS


