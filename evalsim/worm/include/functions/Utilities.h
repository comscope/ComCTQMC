#ifndef EVALSIM_INCLUDE_FUNCTIONS_UTILITIES_H
#define EVALSIM_INCLUDE_FUNCTIONS_UTILITIES_H

#include <tuple>
#include <vector>
#include <string>

#include "Measurements.h"

#include "../../../partition/Functions.h"

#include "../../../../include/measurements/Measurements.h"
#include "../../../../include/measurements/Measurements.h"
#include "../../../../include/io/Matrix.h"
#include "../../../../ctqmc/include/Utilities.h"


namespace evalsim {
    
    namespace worm {
        
        namespace func {
            
            struct iOmega {
                iOmega() = delete;
                iOmega(double const beta, int const is_fermionic=1);
                std::complex<double> operator()(int n) const;
            private:
                double const beta_;
                int const is_fermionic_;
            };
            
        
            struct OmegaMap{
                
                OmegaMap() = delete;
                OmegaMap(std::size_t const n, bool bosonic, bool symmetric);
                
                int operator()(std::size_t const i) const { return val_[i];}
                inline std::size_t pos(int const i) const { auto it = pos_.find(i); return pos_.end() == it ? -1 : it->second;}
                
            private:
                
                int const n_;
                
                std::map<int,std::size_t> pos_;
                std::vector<int> val_;
                
            };
        
        template <typename Value>
        void extend_hyb(std::vector<io::cmat>& hyb, std::vector<io::Matrix<Value>> const& hybMoments, int const nf, double const beta);
        
        template <typename Value>
        void green_and_self_function_on_full_axis(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix,
                                                  std::vector<io::cmat> & hyb, std::vector<io::Matrix<Value>> const& hybMoments, int const nf,
                                                  std::vector<io::cmat>& green, std::vector<io::cmat>& self);
            
        template<typename Value>
        std::vector<io::cmat> get_green_from_obs(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name);
        
        
        template<typename Value>
        std::vector<io::Matrix<Value>> get_moments_from_obs(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name);
        
        }
        
    }
    
}

#endif //EVALSIM










