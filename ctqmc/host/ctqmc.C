

#include <stdexcept>
#include <mpi.h>
#include <memory>
#include <algorithm>


#include "Algebra.h"
#include "../include/MonteCarlo.h"
#include "../../include/parameters/Initialize.h"
#include "../../include/mpi/Utilities.h"


ut::Beta ut::beta;

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    try {
        if(argc != 2) throw std::runtime_error("ctqmc: Wrong number of input parameters!");
        
        std::time_t time;  mpi::cout = mpi::cout_mode::one;
        
        mpi::cout << "Start task at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl << std::endl;
        
        auto const case_name = argv[1];
        
        jsx::value jParams = mpi::read(std::string(std::string(case_name)) + ".json");  params::initialize(jParams); params::complete_worms(jParams);
        if (jParams.is("restart") and jParams("restart").boolean()) jParams["measurements"] = mpi::read(std::string(std::string(case_name))+".meas.json");
        
        jsx::value jSimulation = jsx::array_t{
            jsx::object_t{{ "id", mpi::rank() }, { "config", jsx::read("config_" + std::to_string(mpi::rank()) + ".json", jsx::object_t()) }}
        };
        
        if(jParams.is("complex") ? jParams("complex").boolean() : false) {
            mc::zhMonteCarlo(jParams, jSimulation);
        } else {
            mc::dhMonteCarlo(jParams, jSimulation);
        }
        
        jsx::write(jSimulation("configs")(0), "config_" + std::to_string(mpi::rank()) + ".json");

        mpi::write(jSimulation("measurements"), std::string(std::string(case_name)) + ".meas.json");
        mpi::write(jSimulation("info"),         std::string(std::string(case_name)) + ".info.json");
        
        if(jSimulation.is("error")) mpi::write(jSimulation("error"), std::string(std::string(case_name)) + ".err.json");
        if(jSimulation.is("resample")) jsx::write(jSimulation("resample"), std::string(std::string(case_name)) + ".meas" + std::to_string(mpi::rank()) + ".json");
        
        mpi::cout << "Task of worker finished at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
    }
    catch (std::exception& exc) {
        std::cerr << exc.what() << " ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        return -1;
    }
    catch (...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -2);
#endif
        return -2;
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}












