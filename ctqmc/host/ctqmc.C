#ifndef CTQMC_DRIVER
#define CTQMC_DRIVER

#include "ctqmc.h"

namespace ctqmc{
    //specialized gpu/cpu
    void ctqmc_driver(const char* case_name){
        
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
        if(jSimulation.is("variance")) mpi::write(jSimulation("variance"), std::string(std::string(case_name)) + ".var.json");
        if(jSimulation.is("resample")) jsx::write(jSimulation("resample"), std::string(std::string(case_name)) + ".meas" + std::to_string(mpi::rank()) + ".json");
        
        mpi::cout << "Task of worker finished at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
    }
}

#endif
