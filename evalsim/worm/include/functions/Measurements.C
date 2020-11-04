#include "Measurements.h"

namespace evalsim {
    
    namespace worm {
        
        namespace meas {
            
            namespace impl {
                
                io::cvec read_matsubara_function(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jFunction, bool scale){
                    double const beta = jParams("beta").real64();
                    
                    io::cvec function = jsx::at<io::cvec>(jFunction);
                    
                    if (scale) for (auto& x : function) x*=-1./beta;
                    
                    return function;
                }
                
                io::cvec collapse_antisymmetric_function(io::cvec const& function){
                    
                    io::cvec as_function(function.size()/2);
                    
                    std::size_t start = function.size()/2;
                    
                    auto it_forward = function.begin()+start;
                    auto it_backward = function.begin()+start-1;
                    
                    for (auto& x : as_function){
                        x = 0.5 * ( *it_forward + std::conj(*it_backward) );
                        --it_backward; ++it_forward;
                    }
                    
                    return as_function;
                }
                
            }
        }
    }
}

