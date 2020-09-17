
#include <memory>
#include <algorithm>

#include "../ctqmc/host/Algebra.h"
#include "../ctqmc/device/planar_complex/Algebra.h"

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
        
        mpi::write(jObservables0, std::string(case_name) + ".obs.json");
        
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


extern "C" int CTQMCDriverStart(char* case_name)
{
    try {
        std::time_t time;
        
        mpi::cout = mpi::cout_mode::one;
        
        mpi::cout << "Start task at " << std::ctime(&(time = std::time(nullptr))) << std::endl;
        
        jsx::value jParams = mpi::read(std::string(case_name) + ".json"); params::initialize(jParams); params::complete_worms(jParams);
        std::size_t const streamsPerProcess  = jParams("sim per device").int64();
        std::size_t const processesPerDevice = 1;
        
        std::vector<char> nodeNames = mpi::processor_name();
        mpi::gather(nodeNames, mpi::master);
                
        std::vector<int> deviceId(mpi::number_of_workers(), -1);
        std::vector<std::int64_t> mcIds; //Id of Markov chain
        if(mpi::rank() == mpi::master) {
            
            int deviceCount = 0; cudaErrchk(cudaGetDeviceCount(&deviceCount));
            std::map<std::string,int> rank_on_node;
            
            std::int64_t mcId = 0;
            for(int rank = 0; rank < mpi::number_of_workers(); ++rank) {
                std::string const nodeName(&nodeNames[rank*mpi::processor_name_size()], mpi::processor_name_size());
                
                if(!rank_on_node.count(nodeName)) rank_on_node[nodeName] = 0;
                else rank_on_node[nodeName]++;
                    
                if(rank_on_node[nodeName] < deviceCount)
                    deviceId[rank] = rank_on_node[nodeName];
                
                //mpi::cout << rank_on_node[nodeName] << " " << nodeName << " " << rank << " " << deviceId[rank] << "\n";
                
                mcIds.push_back(mcId);
                mcId += (deviceId[rank] != -1 and streamsPerProcess) ? streamsPerProcess : 1;
                mcIds.push_back(mcId);
            }
        }
        mpi::bcast(deviceId, mpi::master);  mpi::scatter(mcIds, 2, mpi::master);
        
        mpi::cout << mpi::rank() << " " << deviceId[mpi::rank()] << "\n";
        
        jsx::value jSimulation = jsx::array_t(mcIds.back() - mcIds.front());
        for(auto mcId = mcIds.front(); mcId < mcIds.back(); ++mcId)
            jSimulation(mcId - mcIds.front()) = jsx::object_t{
                { "id", mcId },
                { "config", jsx::read("config_" + std::to_string(mcId) + ".json", jsx::object_t()) }
            };
        
        std::int64_t mode = (deviceId[mpi::rank()] != -1 and streamsPerProcess > 0) ? 1 : 0;
        if(mode) {
            //std::cout << "Rank " << mpi::rank() << " gets simulations [" << mcIds.front() << ", " << mcIds.back() << ") and uses device " << deviceId[mpi::rank()] <<  std::endl;
            
            imp::init_device(deviceId[mpi::rank()], processesPerDevice);
            if(jParams.is("complex") ? jParams("complex").boolean() : false) {
                mc::montecarlo<imp::Device, ut::complex>(jParams, jSimulation);
                mc::statistics<ut::complex>(jParams, jSimulation);
            } else {
                mc::montecarlo<imp::Device, double>(jParams, jSimulation);
                mc::statistics<double>(jParams, jSimulation);
            }
            
            
            imp::release_device();
        } else {
            //std::cout << "Rank " << mpi::rank() << " gets simulations [" << mcIds.front() << ", " << mcIds.back() << ") and uses host" << std::endl;
            
            if(jParams.is("complex") ? jParams("complex").boolean() : false) {
                mc::montecarlo<imp::Host, ut::complex>(jParams, jSimulation);
                mc::statistics<ut::complex>(jParams, jSimulation);
            } else {
                mc::montecarlo<imp::Host, double>(jParams, jSimulation);
                mc::statistics<double>(jParams, jSimulation);
            }
        }
        mpi::reduce<mpi::op::sum>(mode, mpi::master);
        jSimulation["info"]["number of GPUs"] = jsx::int64_t(mode);
        
        std::vector<std::size_t> number = { jSimulation("configs").size() };  mpi::gather(number, mpi::master);
        
        mcIds.clear();
        if(mpi::rank() == mpi::master) {
            std::size_t mcId = 0;
            for(int rank = 0; rank < mpi::number_of_workers(); ++rank) {
                mcIds.push_back(mcId);
                mcId += number[rank];
                mcIds.push_back(mcId);
            }
        }
        mpi::scatter(mcIds, 2, mpi::master);
        
        for(std::size_t mcId = mcIds.front(); mcId < mcIds.back(); ++mcId) {
            std::ofstream file("config_" + std::to_string(mcId) + ".json");
            jsx::write(jSimulation("configs")(mcId - mcIds.front()), file);
            file.close();
        }
        
        mpi::write(jSimulation("measurements"), std::string(case_name) + ".meas.json");
        mpi::write(jSimulation("info"),         std::string(case_name) + ".info.json");
        
        
        if(jSimulation.is("error")) mpi::write(jSimulation("error"), std::string(case_name) + ".err.json");
        if(jSimulation.is("resample")) jsx::write(jSimulation("resample"), std::string(case_name) + ".meas" + std::to_string(mpi::rank()) + ".json");
        
        
        mpi::cout << "Task of worker finished at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl;
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



