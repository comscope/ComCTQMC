

#include "Algebra.h"
#include "../include/MonteCarlo.h"

namespace mc{
 
    void zhMonteCarlo(jsx::value const& jParams, jsx::value& jSimulation){
        mc::montecarlo<imp::Host, ut::complex>(jParams, jSimulation);
        mc::statistics<ut::complex>(jParams, jSimulation);
    }
    
}
