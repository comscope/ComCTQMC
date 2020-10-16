#include "libCTQMC.h"

ut::Beta ut::beta;

extern "C" int EvalSim_Main(const char* case_name)
{
    try {
        evalsim::evalsim_driver(case_name);
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

extern "C" int CTQMCDriverStart(const char* case_name)
{
        try {
            ctqmc::ctqmc_driver(case_name);
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

