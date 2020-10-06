
#include <memory>
#include <algorithm>

#include "../ctqmc/host/Algebra.h"

#include "../evalsim/Evalsim.h"

#include "../ctqmc/include/MonteCarlo.h"

#include "../include/parameters/Initialize.h"


ut::Beta ut::beta;



// Interface routine:
// Should be invoked with the case_name, <case_name>.json should exist in the current directory
// It should contains the fields defined by Patrique, including "Atomic" and "QN", although these file will be generated.
//

// adler: this was a main in the original code, we turned it into a named function
//        to be able to contain it in the same library as other mains

jsx::value get_observables(jsx::value const& jParams, std::string const name) {
    jsx::value jMeasurements = mpi::read(name);  io::from_tagged_json(jMeasurements);
    
    if(jParams.is("complex") ? jParams("complex").boolean() : false)
        return evalsim::evalsim<ut::complex>(jParams, jMeasurements);
    else
        return evalsim::evalsim<double>(jParams, jMeasurements);
}

extern "C" int EvalSim_Main(const char* case_name)
{
    try {

        std::time_t time;
        mpi::cout = mpi::cout_mode::one;
        mpi::cout << "Start post-processing at " << std::ctime(&(time = std::time(nullptr))) << std::endl;
        
        jsx::value jParams = mpi::read(std::string(case_name) + ".json"); params::initialize(jParams); params::complete_worms(jParams);
           
        jsx::value jObservables0 = get_observables(jParams, std::string(case_name) + ".meas.json");
           
        std::size_t number_of_mpi_processes = mpi::read(std::string(case_name) + ".info.json")("number of mpi processes").int64();
           
        if(number_of_mpi_processes > 1 && jParams.is("error") && jParams("error").string() == "serial") {
            meas::Error error;
                      
            for(std::size_t i = 0; i < number_of_mpi_processes; ++i)
                    error.add(get_observables(jParams, std::string(case_name) + ".meas" + std::to_string(i) + ".json"), jObservables0);
           
            mpi::write(error.finalize(number_of_mpi_processes, jObservables0), std::string(case_name) + ".err.json");
        }
            
        mpi::write(jObservables0, std::string(std::string(case_name)) + ".obs.json");

        mpi::cout << "End post-processing at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
           
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

    //#ifdef HAVE_MPI
    //    MPI_Init(&argc, &argv);
    //#endif
        try {
            //if(argc != 2) throw std::runtime_error("ctqmc: Wrong number of input parameters!");
            
            std::time_t time;  mpi::cout = mpi::cout_mode::one;
            
            mpi::cout << "Start task at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl << std::endl;
            
            jsx::value jParams = mpi::read(std::string(std::string(case_name)) + ".json");  params::initialize(jParams); params::complete_worms(jParams); 
            if (jParams.is("restart") and jParams("restart").boolean()) jParams["measurements"] = mpi::read(std::string(std::string(case_name))+".meas.json");
            
            jsx::value jSimulation = jsx::array_t{
                jsx::object_t{{ "id", mpi::rank() }, { "config", jsx::read("config_" + std::to_string(mpi::rank()) + ".json", jsx::object_t()) }}
            };
            
            if(jParams.is("complex") ? jParams("complex").boolean() : false) {
                mc::montecarlo<imp::Host, ut::complex>(jParams, jSimulation);
                mc::statistics<ut::complex>(jParams, jSimulation);
            } else {
                mc::montecarlo<imp::Host, double>(jParams, jSimulation);
                mc::statistics<double>(jParams, jSimulation);
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
            
    //#ifdef HAVE_MPI
    //        MPI_Abort(MPI_COMM_WORLD, -1);
    //#endif
            return -1;
        }
        catch (...) {
            std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
            
    //#ifdef HAVE_MPI
    //        MPI_Abort(MPI_COMM_WORLD, -2);
    //#endif
            return -2;
        }
        
    //#ifdef HAVE_MPI
    //    MPI_Finalize();
    //#endif

    return 0;
}



