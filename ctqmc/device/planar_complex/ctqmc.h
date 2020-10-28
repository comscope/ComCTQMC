#ifndef CTQMC_GPU_DRIVER
#define CTQMC_GPU_DRIVER

#include <stdexcept>
#include <mpi.h>

#include "../../host/Algebra.h"
#include "Algebra.h"

#include "../../../evalsim/Evalsim.h"
#include "../../include/MonteCarlo.h"
#include "../../../include/parameters/Initialize.h"

namespace ctqmc{
    void ctqmc_gpu_driver(const char * case_name);
}
#endif













