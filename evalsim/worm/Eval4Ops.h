#ifndef EVALSIM_WORM_EVAL4OPS
#define EVALSIM_WORM_EVAL4OPS


#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

namespace evalsim {
    
    namespace worm {
        
        jsx::value read_4_ops(jsx::value const& jParams, jsx::value const& jMeasurements) {
            double const beta = jParams("beta").real64();
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");

            jsx::value jFunctions;

            for(std::size_t cDagg1 = 0; cDagg1 < 2*jHybMatrix.size(); ++cDagg1)
                for(std::size_t c1 = 0; c1 < 2*jHybMatrix.size(); ++c1)
                    for(std::size_t cDagg2 = 0; cDagg2 < 2*jHybMatrix.size(); ++cDagg2)
                        for(std::size_t c2 = 0; c2 < 2*jHybMatrix.size(); ++c2) {
                            auto name = std::to_string(c1) + "_" + std::to_string(cDagg1) + "_" + std::to_string(cDagg2) + "_" + std::to_string(c2);
                            if(jMeasurements.is(name)) {
                                auto function = jsx::at<io::cvec>(jMeasurements(name));
                                
                                for(auto& entry : function) entry /= beta;
                                
                                jFunctions[name] = std::move(function);
                            };
                        };
            
            return jFunctions;
        };
        
        
        template<typename Value, typename W>
        jsx::value eval_4_ops(ut::wrap<W>, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            
            if(jParams(W::name())("basis").string() == "matsubara") {
                return read_4_ops(jParams, jMeasurements);
            } else if(jParams(W::name())("basis").string() == "legendre") {
                throw std::runtime_error("evalsim::worm::eval_4_ops: legendre not implemented");
            } else
                throw std::runtime_error("evalsim: unknow basis option for green worm");

        }
        
    }
    
}


#endif









