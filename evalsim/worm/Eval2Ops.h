#ifndef EVALSIM_WORM_EVAL2OPS
#define EVALSIM_WORM_EVAL2OPS


#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

namespace evalsim {
    
    namespace worm {
        
        jsx::value read_2_ops(jsx::value const& jParams, jsx::value const& jMeasurements) {
            double const beta = jParams("beta").real64();
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            
            jsx::value jFunctions;
            
            for(std::size_t cDagg = 0; cDagg < 2*jHybMatrix.size(); ++cDagg)
                for(std::size_t c = 0; c < 2*jHybMatrix.size(); ++c) {
                    auto name = std::to_string(c) + "_" + std::to_string(cDagg);
                    if(jMeasurements.is(name)) {
                        auto function = jsx::at<io::cvec>(jMeasurements(name));
                        
                        for(auto& entry : function) entry /= -beta;
                        
                        jFunctions[name] = std::move(function);
                    };
                };
            
            return jFunctions;
        };
        
        
        template<typename Value, typename W>
        jsx::value eval_2_ops(ut::wrap<W>, jsx::value const& jParams, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value& jObservables) {
            
            if(jParams(W::name())("basis").string() == "matsubara") {
                return read_2_ops(jParams, jMeasurements);
            } else if(jParams(W::name())("basis").string() == "legendre") {
                throw std::runtime_error("evalsim::worm::eval_2_ops: legendre not implemented");
            } else
                throw std::runtime_error("evalsim: unknow basis option for green worm");
            
        }
        
    }
    
}


#endif









