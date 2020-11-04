#ifndef CTQMC_INCLUDE_MONTECARLO_H
#define CTQMC_INCLUDE_MONTECARLO_H

#include <ratio>
#include <chrono>
#include <ctime>
#include <tuple>
#include <random>

#include "Params.h"

#include "markovchain/Scheduler.h"
#include "markovchain/MarkovChain.h"

#include "updates/Setup.h"

#include "observables/Observables.h"
#include "observables/Setup.h"

#include "../../include/io/Tag.h"
#include "../../include/measurements/Error.h"
#include "../../evalsim/Evalsim.h"

namespace mc {
    
    template<typename Mode, typename Value>
    void montecarlo(jsx::value jParams, jsx::value& jSimulation);
    
    
    
    template<typename Value>
    void statistics(jsx::value jParams, jsx::value& jSimulation);

}

#include "MonteCarlo.impl.h"

#endif











