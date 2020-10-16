#ifndef EVALSIM_WORM_HEDIN_IMPL
#define EVALSIM_WORM_HEDIN_IMPL

#include "Hedin.h"

namespace evalsim {
    
    namespace worm {
        
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::hedin_ph::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::hedin_ph::Worm::name();
        jsx::value jWorm = jParams(name);
        
        ////Initialization
        func::Frequencies<Value> frequencies(jWorm);
        
        mpi::cout << "Reading hybridisation function ... " << std::flush;
        
        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
        std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
        
        std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
        
        mpi::cout << "Ok" << std::endl;
        
        mpi::cout << "Reading " << name << " function ... " << std::flush;
        
        std::vector<io::ctens> hedin = meas::read_tensor_functions<Value,Fermion,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
        
        mpi::cout << "Ok" << std::endl;
        
        
        mpi::cout << "Gathering one-particle green's function ... " << std::flush;
        
        auto const green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
        auto const green = func::green_function_on_full_axis( green_half_axis );
        func::OmegaMap greenOM(green.size(),false,true);
        
        mpi::cout << "Ok" << std::endl;
        
        
        mpi::cout << "Calculating susceptibility ... " << std::flush;
        
        func::hedin::ph::compute_and_subtract_disconnected<Value>(jParams, frequencies, green, greenOM, jObservables("partition")("occupation"), hedin);
        
        mpi::cout << "Ok" << std::endl;
        
        
        mpi::cout << "Enforcing symmetries ... " << std::flush;
        
        std::vector<io::ctens> susc_symm(hedin.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
        func::hedin::ph::enforce_symmetries<Value>(jParams, frequencies, hedin, susc_symm);
        
        mpi::cout << "Ok" << std::endl;
        
        
        jsx::value jObservablesOut;
        
        jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
        
        return jObservablesOut;
    }
    
    template<typename Value>
    jsx::value evalsim(ut::wrap<cfg::hedin_pp::Worm>, jsx::value jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables) {
        
        auto const name = cfg::hedin_pp::Worm::name();
        jsx::value jWorm = jParams(name);
        
        
        ////Initialization
        func::Frequencies<Value> frequencies(jWorm);
        
        mpi::cout << "Reading hybridisation function ... " << std::flush;
        
        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
        std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
        
        std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
        
        mpi::cout << "Ok" << std::endl;
        
        mpi::cout << "Reading " << name << " function ... " << std::flush;
        
        std::vector<io::ctens> hedin = meas::read_tensor_functions<Value,Fermion,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
        
        mpi::cout << "Ok" << std::endl;
        
        
        mpi::cout << "Gathering one-particle green's function ... " << std::flush;
        
        auto const green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
        auto const green = func::green_function_on_full_axis( green_half_axis );
        func::OmegaMap greenOM(green.size(),false,true);
        
        mpi::cout << "Ok" << std::endl;
        
        
        mpi::cout << "Calculating susceptibility ... " << std::flush;
        
        func::hedin::pp::compute_and_subtract_disconnected<Value>(jParams, frequencies, green, greenOM, jObservables("partition")("occupation"), hedin);
        
        mpi::cout << "Ok" << std::endl;
        
        
        mpi::cout << "Enforcing symmetries ... " << std::flush;
        
        std::vector<io::ctens> susc_symm(hedin.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
        func::hedin::pp::enforce_symmetries<Value>(jParams, frequencies, hedin, susc_symm);
        
        mpi::cout << "Ok" << std::endl;
        
        
        jsx::value jObservablesOut;
        
        jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
        
        return jObservablesOut;
    }
        
        
    }
    
}


#endif









