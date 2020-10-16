#ifndef CTQMC_INCLUDE_BATH_BATH_IMPL_H
#define CTQMC_INCLUDE_BATH_BATH_IMPL_H

#include "Bath.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace bath {
    
template<typename Value>
template<typename Type>
void Bath<Value>::add(Hyb<Value> const& hyb, Type type) {
    if(update_.get() == nullptr)
        update_.reset(new Update<Value, Type>());
    else if(update_->type().val_ != get_type_id<Type>().val_)
        throw std::runtime_error("bath::add: update type missmatch");
    
    static_cast<Update<Value, Type>&>(*update_).add(type, *this, hyb);
};

}

#endif
