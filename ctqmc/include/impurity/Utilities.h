#ifndef CTQMC_INCLUDE_IMPURITY_UTILITIES_H
#define CTQMC_INCLUDE_IMPURITY_UTILITIES_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>

#include "../../../include/linalg/LinAlg.h"
#include "../../../include/linalg/Operators.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    template <typename Value>
    io::rvec blockNorms(jsx::value const& jOperator);

    struct NormCollection{
        
        NormCollection() = delete;
        NormCollection(std::size_t const N) : norms_(N), normsDagg_(N){}
        NormCollection(std::vector<io::rvec> norms, std::vector<io::rvec> normsDagg) : norms_(norms), normsDagg_(normsDagg){}
        
        std::vector<io::rvec> const& norms() const {return norms_;}
        std::vector<io::rvec> const& normsDagg() const {return normsDagg_;}
        
    private:
        std::vector<io::rvec> norms_;
        std::vector<io::rvec> normsDagg_;
        
    };

    template <typename Value>
    NormCollection gatherNorms(jsx::value const& jMPI, jsx::value const& jOperators);

}

#include "Utilities.impl.h"

#endif  
