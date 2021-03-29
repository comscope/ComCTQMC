#include "Evalsim.h"

namespace evalsim {


    
    void evalsim_driver(const char* case_name){
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
    
    
    //---------------------------------------------------------------------------------------------------------------------------------------
    
    template <typename Value>
    void complete_params(jsx::value & jParams){
        bool const ising = (jParams("hloc")("two body").is<jsx::object_t>() && jParams("hloc")("two body").is("approximation")) ? (jParams("hloc")("two body")("approximation").string() == "ising") : false;
        
        jParams["hloc"] = ga::read_hloc<Value>("hloc.json");
        
        jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
        
        if(jParams.is("dyn"))
            jParams("dyn")("functions") = mpi::read(jParams("dyn")("functions").string());
        
        jsx::value jPartition = jParams(cfg::partition::Worm::name());
        
        //In case we need the raw input rather than, e.g., the sparse representation
        jParams["partition record"] = jPartition;
        
        opt::complete_qn<Value>(jParams, jPartition["quantum numbers"]);
        
        opt::complete_observables<Value>(jParams, jPartition["observables"], ising);
        
        jParams(cfg::partition::Worm::name()) = jPartition;
        
    }

    jsx::value get_observables(jsx::value & jParams, std::string const name) {
        jsx::value jMeasurements = mpi::read(name);  io::from_tagged_json(jMeasurements);
        
        if(jParams.is("complex") ? jParams("complex").boolean() : false)
            return evalsim::evalsim<ut::complex>(jParams, jMeasurements);
        else
            return evalsim::evalsim<double>(jParams, jMeasurements);
    }

    template<typename Value>
    bool add_dynamic(jsx::value& jStatic, jsx::value const& jDynamic)
    {
        if(jStatic.is<io::Vector<Value>>() && jDynamic.is<io::Vector<Value>>()) {
            auto& dyn = jsx::at<io::Vector<Value>>(jDynamic);
            auto& stc = jsx::at<io::Vector<Value>>(jStatic);
            
            if(stc.size() != dyn.size()) return false;
            
            for(std::size_t i = 0; i < dyn.size(); ++i) stc[i] += dyn[i];
            
            return true;
        }
        return false;
    }
    
    template bool add_dynamic<double>(jsx::value& jStatic, jsx::value const& jDynamic);
    template bool add_dynamic<ut::complex>(jsx::value& jStatic, jsx::value const& jDynamic);
    
    void add_dynamics(jsx::value const& jParams, jsx::value& jMeasurements, std::string const worm, std::string const meas)
    {
        if(jParams.is("dyn") && jParams.is(worm + meas)) {
            auto& jDynamic = jMeasurements(worm)("dynamic");
            auto& jStatic = jMeasurements(worm + meas)("static");

            std::set<std::string> entries;
        
            for(auto entry : jDynamic.object()) entries.insert(entry.first);
            for(auto entry : jStatic.object())  entries.insert(entry.first);
        
            for(auto entry : entries) {
                if(jDynamic.is(entry) && !jStatic.is(entry))
                    jStatic[entry] = jDynamic(entry);
                else if (jDynamic.is(entry) && jStatic.is(entry)) {
                    if(!add_dynamic<double>(jStatic(entry), jDynamic(entry)) && !add_dynamic<ut::complex>(jStatic(entry), jDynamic(entry)))
                        throw std::runtime_error("evalsim::add_dynamics: missmatch for " + worm + meas);
                }
            }
        }
    }
    
    template<typename Value>
    jsx::value evalsim(jsx::value & jParams, jsx::value& jMeasurements){
        jsx::value jObservables;
        
        complete_params<Value>(jParams);
        
        mpi::cout << "Begin evaluating partition measurements" << std::endl;
        
        jObservables[cfg::partition::Worm::name()] = partition::evalsim<Value>(jParams, jMeasurements(cfg::partition::Worm::name()));
        
        mpi::cout << "End evaluating partition measurements" << std::endl;
        
        // This option is set in CTQMC if `all errors' is set. The intent is to enable fast error computation.
        // While not all worms are slow to evaluate, they require (in some cases) the occupations.
        // The occupation can be a difficult observable to compute in problems with large invariant subspaces.
        bool const lpp = jParams.is("limited post-processing") ? jParams("limited post-processing").boolean() : false;
        if (!lpp){
            
            add_dynamics(jParams, jMeasurements, cfg::green::name, " impr");
            add_dynamics(jParams, jMeasurements, cfg::green::name, " imprsum");
                
            add_dynamics(jParams, jMeasurements, cfg::vertex::name, " impr");
            add_dynamics(jParams, jMeasurements, cfg::vertex::name, " imprsum");
                
            add_dynamics(jParams, jMeasurements, cfg::hedin_ph::name, " impr");
            add_dynamics(jParams, jMeasurements, cfg::hedin_ph::name, " imprsum");
                
            add_dynamics(jParams, jMeasurements, cfg::hedin_pp::name, " impr");
            add_dynamics(jParams, jMeasurements, cfg::hedin_pp::name, " imprsum");
                
            cfg::for_each_type<cfg::Worm>::apply(worm_clean_functor<Value>(), jParams, jMeasurements);
                
            cfg::for_each_type<cfg::Worm>::apply(worm_evalsim_functor<Value>(), jParams, jMeasurements, jObservables);
                
            if (jParams.is("kernels")){
                jObservables["kernels"] = worm::evaluateKernels<Value>(jParams,jObservables);
                
                if (jParams("kernels").is("full") ? jParams("kernels")("full").boolean() : false)
                    jObservables["Asymptotic Full Vertex"] = worm::evaluateFullVertexFromKernels<Value>(jParams,jObservables);
            }
            
        }
        
        return jObservables;
        
    }
    
    template jsx::value evalsim<double>(jsx::value & jParams, jsx::value& jMeasurements);
    template jsx::value evalsim<ut::complex>(jsx::value & jParams, jsx::value& jMeasurements);
    
}












