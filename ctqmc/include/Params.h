#ifndef CTQMC_INCLUDE_PARAMS_H
#define CTQMC_INCLUDE_PARAMS_H

#include "Utilities.h"
#include "config/Worms.h"

#include "../../include/mpi/Utilities.h"
#include "../../include/atomic/Generate.h"
#include "../../include/options/Options.h"

namespace params {

    template<typename Value>
    void complete_impurity(jsx::value& jParams)
    {
        opt::complete_hloc<Value>(jParams);
        jParams("hloc") = ga::construct_hloc<Value>(jParams("hloc"));
        mpi::write(jParams("hloc"), "hloc.json");
        
        jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
        
        jParams("hybridisation")("functions") = mpi::read(jParams("hybridisation")("functions").string());
        
        if(jParams.is("dyn"))
            jParams("dyn")("functions") = mpi::read(jParams("dyn")("functions").string());
        
        jParams["mpi structure"] = mpi::mpi_structure();
    };
    
    
void complete_worm(jsx::value& jParams, std::string const worm);
    
void complete_worms(jsx::value& jParams);

}

#endif
