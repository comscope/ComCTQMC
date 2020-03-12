#include "Evalsim.h"


jsx::value get_observables(jsx::value const& jParams, std::string const name) {
    jsx::value jMeasurements = mpi::read(name);  io::from_tagged_json(jMeasurements);
    
    if(jParams.is("complex") ? jParams("complex").boolean() : false)
        return evalsim::evalsim<ut::complex>(jParams, jMeasurements);
    else
        return evalsim::evalsim<double>(jParams, jMeasurements);
}


int main(int argc, char** argv)
{
    try {

        if(argc != 2)
            throw std::runtime_error("EvalSim: Wrong number of input parameters !");

        jsx::value jParams = mpi::read(std::string(argv[1]) + ".json");  params::complete_worms(jParams);
        
        
        jsx::value jObservables0 = get_observables(jParams, std::string(argv[1]) + ".meas.json");
        
        
        std::size_t number_of_mpi_processes = mpi::read(std::string(argv[1]) + ".info.json")("number of mpi processes").int64();
        
        if(number_of_mpi_processes > 1 && jParams.is("error") && jParams("error").string() == "serial") {
            meas::Error error;
            
            for(std::size_t i = 0; i < number_of_mpi_processes; ++i)
                error.add(get_observables(jParams, std::string(argv[1]) + ".meas" + std::to_string(i) + ".json"), jObservables0);
 
            mpi::write(error.finalize(number_of_mpi_processes, jObservables0), std::string(argv[1]) + ".err.json");
        }
        
        mpi::write(jObservables0, std::string(argv[1]) + ".obs.json");
        
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












