#ifndef LIBCTQMC
#define LIBCTQMC

#include "../ctqmc/host/ctqmc.h"
#include "../evalsim/Evalsim.h"

extern "C" int CTQMCDriverStart(const char* case_name);
extern "C" int EvalSim_Main(const char* case_name);

#endif
