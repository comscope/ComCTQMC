#ifndef INCLUDE_LINALG_LINALG_IMPL_H
#define INCLUDE_LINALG_LINALG_IMPL_H

#include "LinAlg.h"

namespace linalg {
    
template <typename T>
typename std::enable_if<!std::is_integral<T>::value, std::vector<T>>::type linspace(std::size_t const size, T const begin, T const end, bool endpoint){
    
    std::vector<T> v(size);
    
    if (size){
        
        auto const delta = (end - begin) / static_cast<double>(size - endpoint);
        
        v[0] = begin;
        for (std::size_t i = 1; i<size-endpoint; i++)
        v[i] = v[i-1] + delta;
        
        if (endpoint) v[size-1] = end;
        
    }
    
    return v;
    
}
}

#endif
