#ifndef INCLUDE_MPI_UTILITIES_IMPL_H
#define INCLUDE_MPI_UTILITIES_IMPL_H

#include "Utilities.h"

namespace mpi {		
		
    
    template<typename T>
    std::ostream& operator<<(Cout const& c, T t) {
        if(c.mode_ == cout_mode::every) return std::cout << t;
        if(c.mode_ == cout_mode::one) return rank() == master ? std::cout << t : Cout::null_;
        return Cout::null_;
    };

}


#endif
