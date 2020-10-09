

#include "Algebra.h"
#include "../include/MonteCarlo.h"

namespace mc{
 
    void dhMonteCarlo(jsx::value const& jParams, jsx::value& jSimulation){
        mc::montecarlo<imp::Host, double>(jParams, jSimulation);
        mc::statistics<double>(jParams, jSimulation);
    }
}

