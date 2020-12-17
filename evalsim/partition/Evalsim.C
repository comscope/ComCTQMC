#include "Evalsim.h"

namespace evalsim {
    
    namespace partition {
        
        template<typename Value>
        jsx::value evalsim(jsx::value jParams, jsx::value const& jMeasurements) {

            //This options allows one to not compute those objects which can take a while to compute
            bool const lpp = jParams.is("limited post-processing") ? jParams("limited post-processing").boolean() : false;
            
            //TODO: Move this?
            bool const ising = (jParams("hloc")("two body").is<jsx::object_t>() && jParams("hloc")("two body").is("approximation")) ? (jParams("hloc")("two body")("approximation").string() == "ising") : false;
            
            jParams["hloc"] = ga::read_hloc<Value>("hloc.json");
            
            jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
            
            if(jParams.is("dyn"))
                jParams("dyn")("functions") = mpi::read(jParams("dyn")("functions").string());
            
            
            jsx::value jPartition = jParams(cfg::partition::Worm::name());
            
            opt::complete_qn<Value>(jParams, jPartition["quantum numbers"]);
            
            opt::complete_observables<Value>(jParams, jPartition["observables"],ising);
            
            jsx::value jObservables;
            
            
            
            jObservables["sign"] = jMeasurements("sign");
            
            jObservables["expansion histogram"] = jMeasurements("expansion histogram");
            
            jObservables["scalar"] = get_scalar<Value>(jParams, jPartition, jMeasurements);
            

            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            io::Matrix<Value> occupation(jHybMatrix.size(), jHybMatrix.size()), correlation(jHybMatrix.size(), jHybMatrix.size());
            if (!lpp)
                jObservables["occupation"] = get_occupation<Value>(jParams, jMeasurements, occupation, correlation);
            
            
            
            mpi::cout << "Reading hybridisation function ... " << std::flush;
            
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = func::get_hybridisation<Value>(jParams);
            
            mpi::cout << "Ok" << std::endl;
            
            
            
            mpi::cout << "Reading green function ... " << std::flush;
            
            std::vector<io::cmat> green = meas::read_functions<Value>(jMeasurements("green"), jParams, jPartition, hyb.size());
            
            mpi::cout << "Ok" << std::endl;
            
            
            
            mpi::cout << "Calculating self-energy with dyson ... " << std::flush;
            
            std::vector<io::cmat> selfDyson = func::get_self_dyson<Value>(jParams, green, hyb);
            
            mpi::cout << "OK" << std::endl;
            
            
            
            std::vector<io::cmat> selfenergy;
            
            if(jPartition.is("green bulla") ? jPartition("green bulla").boolean() : true) {
                
                mpi::cout << "Calculating self-energy with bulla ... " << std::flush;
                
                selfenergy.resize(green.size(), io::cmat(jParams("hybridisation")("matrix").size(), jParams("hybridisation")("matrix").size()));
                
                auto bulla = meas::read_conj_functions<Value>(jMeasurements("bullaL"), jMeasurements("bullaR"), jParams, jPartition, hyb.size());
                
                for(std::size_t n = 0; n < green.size(); ++n) {
                    io::cmat green_inv = linalg::inv(green[n]);
                    
                    linalg::mult<ut::complex>('n', 'n', .5, green_inv, bulla.first[n], .0, selfenergy[n]);
                    linalg::mult<ut::complex>('n', 'n', .5, bulla.second[n], green_inv, 1., selfenergy[n]);
                }
                
                mpi::cout << "Ok" << std::endl;

            } else
                selfenergy = std::move(selfDyson);
            
            
            if (!lpp){
                mpi::cout << "Calculating green moments ... " << std::flush;
                
                std::vector<io::Matrix<Value>> greenMoments = get_green_moments(jParams, hybMoments, jMeasurements, jObservables("scalar"));
                
                mpi::cout << "OK" << std::endl;
                
                
                
                mpi::cout << "Calculating self-energy moments ... " << std::flush;
                
                std::vector<io::Matrix<Value>> selfMoments = get_self_moments(jParams, hybMoments, greenMoments);
                
                mpi::cout << "OK" << std::endl;
            
            
                mpi::cout << "Adding self-energy high frequency tail ... "  << std::flush;
                
                auto tail_length = hyb.size();
                if (jParams.is("analytical continuation"))
                    tail_length = std::max(tail_length, static_cast<std::size_t>(1.5*jParams("analytical continuation")("nf").int64()));
                
                func::add_self_tail(jParams, selfenergy, selfMoments, tail_length);   //scheisse Value !! Allgemein scheiss moments ...
                    
                jObservables["self-energy"] = func::write_functions(jParams, selfenergy, selfMoments);
                    
                if(selfDyson.size()) jObservables["self-energy-dyson"] = func::write_functions(jParams, selfDyson, selfMoments);

                mpi::cout << "Ok" << std::endl;
                
                
                mpi::cout << "Adding green function high frequency tail ... " << std::flush;
                
                func::add_green_tail<Value>(jParams, hyb, selfenergy, green);
                
                jObservables["green"] = func::write_functions(jParams, green, greenMoments);
                
                mpi::cout << "Ok" << std::endl;
                
                
                if(jParams.is("analytical continuation")){
                    
                    auto const aux = func::get_aux_green(jParams, selfenergy, selfMoments);
                    
                    auto const jAux = func::write_functions<Value>(jParams, aux);
                    jObservables["aux green matsubara"] = jAux;
                    
                    std::size_t const nf = jParams("analytical continuation")("nf").int64();
                    
                    std::size_t ntau = jParams("analytical continuation").is("ntau") ?
                        jParams("analytical continuation")("ntau").int64() :
                        static_cast<std::size_t>(mpi::number_of_workers()/2);
                    
                    jObservables["aux green"] = func::fourier_transform(jParams, jAux, tail_length, nf, ntau);
                    
                }
                
            } else {
                
                jObservables["self-energy"] = func::write_functions<Value>(jParams, selfenergy);
                    
                jObservables["green"] = func::write_functions<Value>(jParams, green);
                
            }
            
            
            if(jPartition.is("quantum number susceptibility") ? jPartition("quantum number susceptibility").boolean() : false)
                
                jObservables["susceptibility"] = get_qn_susc(jParams, jPartition, jMeasurements, jObservables("scalar"));
            
            if(!lpp)
                if((jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false) ||
                   (jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false)) {
                    
                    mpi::cout << "Calculating occupation susceptibility moments ... " << std::flush;
                    
                    io::rmat moments = get_occupation_susc_moments<Value>(jParams, jPartition, jMeasurements, jObservables);
                    
                    mpi::cout << "Ok" << std::endl;
                    
                    
                    if(jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false) {
                        
                        mpi::cout << "Reading bulla occupation susceptibility ... " << std::flush;
                        
                        jObservables["occupation-susceptibility-bulla"] = get_occupation_susc_bulla(jParams, jPartition, jMeasurements, moments, occupation, correlation, jObservables);
                    
                        mpi::cout << "Ok" << std::endl;
                        
                    }
                    
                    
                    if(jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false) {
                        
                        mpi::cout << "Reading direct occupation susceptibility ... " << std::flush;
                        
                        jObservables["occupation-susceptibility-direct"] = get_occupation_susc_direct(jParams, jPartition, jMeasurements, moments, occupation, correlation, jObservables);
                        
                        mpi::cout << "Ok" << std::endl;
                        
                    }
        
                }

            if(jPartition.is("probabilities"))
                
                jObservables["probabilities"] = get_probabilities<Value>(jParams, jPartition, jMeasurements);
            
            
            if(jPartition.is("print eigenstates") and jPartition("print eigenstates").boolean())
                
                jObservables["eigenstates"] = get_eigenstates<Value>(jParams, jPartition, jMeasurements);
            
            
            
            if(jPartition.is("print density matrix") and jPartition("print density matrix").boolean())
                
                jObservables["density matrix"] = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            
            
            
            return jObservables;
        }
        
        template jsx::value evalsim<double>(jsx::value jParams, jsx::value const& jMeasurements);
        template jsx::value evalsim<ut::complex>(jsx::value jParams, jsx::value const& jMeasurements);
        
    }
    
}

    
