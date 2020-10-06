#ifndef EVALSIM_DRIVER
#define EVALSIM_DRIVER

#include <memory>
#include <algorithm>

#include "../ctqmc/host/Algebra.h"
#include "../evalsim/Evalsim.h"
#include "../ctqmc/include/MonteCarlo.h"
#include "../include/parameters/Initialize.h"

extern "C" int EvalSim_Main(const char* case_name);

jsx::value get_observables(jsx::value const& jParams, std::string const name) {
    jsx::value jMeasurements = mpi::read(name);  io::from_tagged_json(jMeasurements);
    
    if(jParams.is("complex") ? jParams("complex").boolean() : false)
        return evalsim::evalsim<ut::complex>(jParams, jMeasurements);
    else
        return evalsim::evalsim<double>(jParams, jMeasurements);
}

//not specialized gpu/cpu
void EVALSIM(const char* case_name){
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

#endif
