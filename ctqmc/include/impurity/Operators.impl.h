#ifndef CTQMC_INCLUDE_IMPURITY_OPERATORS_IMPL_H
#define CTQMC_INCLUDE_IMPURITY_OPERATORS_IMPL_H

#include "Operators.h"

namespace imp {
    
    template<typename Mode, typename Value>
    template<typename... Args>
    Matrix<Mode, Value>& Operator<Mode,Value>::mat(int s, Args&& ... args) {
        if(isMat(s)) throw std::runtime_error("imp::Operator::mat: matrix allocated");
        new(mat_ + s) Matrix<Mode, Value>(std::forward<Args>(args)...);
        isMat_.set(s); return mat_[s];
    }

}

#endif  
