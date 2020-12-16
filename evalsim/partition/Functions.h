#ifndef EVALSIM_PARTITION_FUNCTIONS
#define EVALSIM_PARTITION_FUNCTIONS

#include <tuple>

#include "../../include/measurements/Measurements.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/atomic/Generate.h"
#include "../../ctqmc/include/bath/Hyb.h"
#include "../../ctqmc/include/Utilities.h"


namespace evalsim {
    
    namespace partition {
        
        namespace func {
            
            struct iOmega {
                iOmega() = delete;
                iOmega(double beta);
                std::complex<double> operator()(int n) const;
            private:
                double const beta_;
            };
        
            
            template<typename Value>
            io::Matrix<Value> get_matrix(std::map<std::string, Value> const& entries, jsx::value const& jMatrix);
            
            
            template<typename Value>
            std::vector<io::Matrix<Value>> get_function_matrix(std::map<std::string, io::Vector<Value>> const& functions, jsx::value const& jMatrix);
            
            
            template<typename Value>
            std::map<std::string, Value> get_entries(io::Matrix<Value> const& matrix, jsx::value const& jMatrix);
            
            
            template<typename Value>
            std::map<std::string, io::Vector<Value>> get_function_entries(std::vector<io::Matrix<Value>> const& functionMatrix, jsx::value const& jMatrix);
            

            template<typename Value>
            std::vector<io::Matrix<Value>> get_hybridisation_moment(std::map<std::string, io::cvec> const& functions, jsx::value const& jParams, jsx::value const& jHybMatrix);
            
            
            template<typename Value>
            std::tuple<std::vector<io::cmat>, std::vector<io::Matrix<Value>>> get_hybridisation(jsx::value const& jParams);
            
            
            template<typename Value>
            std::vector<io::cmat> get_self_dyson(jsx::value const& jParams, std::vector<io::cmat> const& green, std::vector<io::cmat> const& hyb);
            
            
            template<typename Value>
            void add_self_tail(jsx::value const& jParams, std::vector<io::cmat>& function, std::vector<io::Matrix<Value>> const& moments, std::size_t hybSize);
            
            
            template<typename Value>
            void add_green_tail(jsx::value const& jParams, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& selfenergy, std::vector<io::cmat>& green);
            
            
            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix, std::vector<io::Matrix<Value>> const& momentsMatrix);
        
            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix);
            
            template<typename Value>
            std::vector<io::cmat> get_aux_green(jsx::value const& jParams, std::vector<io::cmat> const& selfenergy, std::vector<io::Matrix<Value>> const& selfMoments, std::vector<io::cmat> const& hyb);
                
        }
        
    }
    
}

#endif //EVALSIM










