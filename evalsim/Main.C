#include "Evalsim.h"

ut::Beta ut::beta;

int main(int argc, char** argv)
{
    #ifdef HAVE_MPI
        MPI_Init(&argc, &argv);
    #endif
    try {

        if(argc != 2)
            throw std::runtime_error("EvalSim: Wrong number of input parameters !");
        
        evalsim::evalsim_driver(argv[1]);
    }
    catch(std::exception& exc) {
        std::cerr << exc.what() << "( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        
        return -1;
    }
    catch(...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        
        return -2;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return 0;
}












