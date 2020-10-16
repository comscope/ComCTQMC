#ifndef CTQMC_INCLUDE_PARAMS_H
#define CTQMC_INCLUDE_PARAMS_H

#include "Utilities.h"
#include "config/Worms.h"

#include "../../include/mpi/Utilities.h"
#include "../../include/atomic/Generate.h"
#include "../../include/options/Options.h"

namespace params {

    template<typename Value>
void complete_impurity(jsx::value& jParams);
    
    
void complete_worm(jsx::value& jParams, std::string const worm);
    
void complete_worms(jsx::value& jParams);

}

#include "Params.impl.h"

#endif
