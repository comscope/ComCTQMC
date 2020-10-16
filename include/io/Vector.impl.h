#ifndef INCLUDE_IO_VECTOR_IMPL_H
#define INCLUDE_IO_VECTOR_IMPL_H

#include "Vector.h"

namespace io {
  
template<typename T, typename std::enable_if<!std::is_same<std::complex<double>, T>::value, int>::type>   // sfinae so that it is taken for ivec and rvec. However, not sure anymore why function overload below is not sufficient
jsx::value encode(std::vector<T> const& source, bool b64) {
    if(b64) return base64::encode(source);
    
    return jsx::array_t(source.begin(), source.end());
};

template<typename T, typename std::enable_if<!std::is_same<std::complex<double>, T>::value, int>::type>
void decode(jsx::value const& source, std::vector<T>& dest) {
    if(source.is<jsx::array_t>()) {
        for(auto const& x : source.array()) dest.push_back(x.real64());
    } else if(source.is<jsx::string_t>()) {
        base64::decode(source.string(), dest);
    } else
        throw std::runtime_error("io::decode: invalid format");
};

};

#endif
