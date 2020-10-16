#ifndef CTQMC_HOST
#define CTQMC_HOST

#include <memory>
#include <algorithm>

#include "Algebra.h"

#include "../../evalsim/Evalsim.h"
#include "../include/MonteCarlo.h"
#include "../../include/parameters/Initialize.h"

namespace ctqmc{

    void ctqmc_driver(const char* case_name);

}

#endif
