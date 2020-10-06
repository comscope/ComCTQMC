
#include <memory>
#include <algorithm>

#include "libCTQMC_GPU.h"
#include "libEVALSIM.h"

ut::Beta ut::beta;

// Interface routine:
// Should be invoked with the case_name, <case_name>.json should exist in the current directory
// It should contains the fields defined by Patrique, including "Atomic" and "QN", although these file will be generated.
//

// adler: this was a main in the original code, we turned it into a named function
//        to be able to contain it in the same library as other mains


extern "C" int EvalSim_Main(const char* case_name)
{
    try {
        EVALSIM(case_name);
    }
    catch(std::exception& exc) {
        std::cerr << exc.what() << "( Thrown from worker " << mpi::rank() << " )" << std::endl;
           
        return -1;
    }
    catch(...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
           
        return -2;
    }
       
    return 0;
}

extern "C" int CTQMCDriverStart(char* case_name)
{
    try {
        CTQMC(case_name);
    }
    catch (ut::out_of_memory error) {
        std::cerr << "Out of memory" << std::endl;
         
        return -1;
    }
    catch (std::exception& exc) {
        std::cerr << exc.what() << " ( Thrown from worker " << mpi::rank() << " )" << std::endl;
          
        return -1;
    }
    catch (...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
            
        return -2;
    }

    return 0;
}



