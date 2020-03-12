#ifndef EVALSIM_PARTITION_EVALSIM
#define EVALSIM_PARTITION_EVALSIM

#include <memory>
#include <algorithm>
#include <complex>

#include "ReadFunctions.h"
#include "Functions.h"

#include "Moments.h"
#include "Scalar.h"
#include "Occupation.h"
#include "Probabilities.h"
#include "QuantumNumberSusc.h"
#include "OccupationSusc.h"

#include "../../include/atomic/Generate.h"
#include "../../include/linalg/Operators.h"
#include "../../include/options/Options.h"

#include "../../ctqmc/include/config/Worms.h"

namespace evalsim {
    
    namespace partition {
        
        template<typename Value>
        jsx::value evalsim(jsx::value jParams, jsx::value const& jMeasurements) {
            
            jParams["hloc"] = ga::read_hloc<Value>("hloc.json");
            
            jParams["operators"] = ga::construct_annihilation_operators<Value>(jParams("hloc"));
            
            if(jParams.is("dyn"))
                jParams("dyn")("functions") = mpi::read(jParams("dyn")("functions").string());
            
            
            
            jsx::value jPartition = jParams(cfg::partition::Worm::name());
            
            opt::complete_qn<Value>(jParams, jPartition["quantum numbers"]);
            
            opt::complete_observables<Value>(jParams, jPartition["observables"]);
            
            for(auto& jObs : jPartition("observables").object())
                jObs.second = ga::construct_observable<Value>(jParams("hloc"), jObs.second);

            
            jsx::value jObservables;
            
            
            
            jObservables["sign"] = jMeasurements("sign");
            
            jObservables["expansion histogram"] = jMeasurements("expansion histogram");
            
            jObservables["scalar"] = get_scalar<Value>(jParams, jPartition, jMeasurements);
            
            
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            io::Matrix<Value> occupation(jHybMatrix.size(), jHybMatrix.size()), correlation(jHybMatrix.size(), jHybMatrix.size());
            
            jObservables["occupation"] = get_occupation<Value>(jParams, jMeasurements, occupation, correlation);
            
            
            
            std::cout << "Reading hybridisation function ... " << std::flush;
            
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            
            std::tie(hyb, hybMoments) = func::get_hybridisation<Value>(jParams);
            
            std::cout << "Ok" << std::endl;
            
            
            
            std::cout << "Reading green function ... " << std::flush;
            
            std::vector<io::cmat> green = meas::read_functions<Value>(jMeasurements("green"), jParams, jPartition, hyb.size());
            
            std::cout << "Ok" << std::endl;
            
            
            
            std::cout << "Calculating self-energy with dyson ... " << std::flush;
            
            std::vector<io::cmat> selfDyson = func::get_self_dyson<Value>(jParams, green, hyb);
            
            std::cout << "OK" << std::endl;
            
            
            
            std::vector<io::cmat> selfenergy;
            
            if(jPartition.is("green bulla") ? jPartition("green bulla").boolean() : true) {
                
                std::cout << "Calculating self-energy with bulla ... " << std::flush;
                
                selfenergy.resize(green.size(), io::cmat(jParams("hybridisation")("matrix").size(), jParams("hybridisation")("matrix").size()));
                
                auto bulla = meas::read_conj_functions<Value>(jMeasurements("bullaL"), jMeasurements("bullaR"), jParams, jPartition, hyb.size());
                
                for(std::size_t n = 0; n < green.size(); ++n) {
                    io::cmat green_inv = linalg::inv(green[n]);
                    
                    linalg::mult<ut::complex>('n', 'n', .5, green_inv, bulla.first[n], .0, selfenergy[n]);
                    linalg::mult<ut::complex>('n', 'n', .5, bulla.second[n], green_inv, 1., selfenergy[n]);
                }
                
                std::cout << "Ok" << std::endl;

            } else
                selfenergy = std::move(selfDyson);
            
            
            
            std::cout << "Calculating green moments ... " << std::flush;
            
            std::vector<io::Matrix<Value>> greenMoments = get_green_moments(jParams, hybMoments, jMeasurements, jObservables("scalar"));
            
            std::cout << "OK" << std::endl;
            
            
            
            std::cout << "Calculating self-energy moments ... " << std::flush;
            
            std::vector<io::Matrix<Value>> selfMoments = get_self_moments(jParams, hybMoments, greenMoments);
            
            std::cout << "OK" << std::endl;
            
            
            
            std::cout << "Adding self-energy high frequency tail ... "  << std::flush;
            
            func::add_self_tail(jParams, selfenergy, selfMoments, hyb.size());   //scheisse Value !! Allgemein scheiss moments ...
                
            jObservables["self-energy"] = func::write_functions(jParams, selfenergy, selfMoments);
                
            if(selfDyson.size()) jObservables["self-energy-dyson"] = func::write_functions(jParams, selfDyson, selfMoments);

            std::cout << "Ok" << std::endl;
            
            
            
            std::cout << "Adding green function high frequency tail ... " << std::flush;
            
            func::add_green_tail<Value>(jParams, hyb, selfenergy, green);
            
            jObservables["green"] = func::write_functions(jParams, green, greenMoments);
            
            std::cout << "Ok" << std::endl;
            

            
            if(jPartition.is("quantum number susceptibility") ? jPartition("quantum number susceptibility").boolean() : false)
                
                jObservables["susceptibility"] = get_qn_susc(jParams, jPartition, jMeasurements, jObservables("scalar"));
            
            
            
            if((jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false) ||
               (jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false)) {
                
                std::cout << "Calculating occupation susceptibility moments ... " << std::flush;
                
                io::rmat moments = get_occupation_susc_moments<Value>(jParams, jPartition, jMeasurements, jObservables);
                
                std::cout << "Ok" << std::endl;
                
                
                if(jPartition.is("occupation susceptibility bulla")  ? jPartition("occupation susceptibility bulla").boolean()  : false) {
                    
                    std::cout << "Reading bulla occupation susceptibility ... " << std::flush;
                    
                    jObservables["occupation-susceptibility-bulla"] = get_occupation_susc_bulla(jParams, jPartition, jMeasurements, moments, occupation, correlation, jObservables);
                
                    std::cout << "Ok" << std::endl;
                    
                }
                
                
                if(jPartition.is("occupation susceptibility direct") ? jPartition("occupation susceptibility direct").boolean() : false) {
                    
                    std::cout << "Reading direct occupation susceptibility ... " << std::flush;
                    
                    jObservables["occupation-susceptibility-direct"] = get_occupation_susc_direct(jParams, jPartition, jMeasurements, moments, occupation, correlation, jObservables);
                    
                    std::cout << "Ok" << std::endl;
                    
                }
    
            }

            
            if(jPartition.is("probabilities"))
                
                jObservables["probabilities"] = get_probabilities<Value>(jParams, jPartition, jMeasurements);
            
            
            
            return jObservables;
        }
        
    }
    
}

    
#endif
    
    
    
    
    
    
    
    
    
