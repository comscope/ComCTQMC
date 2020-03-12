#ifndef EVALSIM_PARTITION_OCCUPATION_H
#define EVALSIM_PARTITION_OCCUPATION_H


#include "ReadDensityMatrix.h"

#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"

namespace evalsim {
    
    namespace partition {
        
        
        template<typename Value>
        jsx::value get_occupation(jsx::value const& jParams, jsx::value const& jMeasurements, io::Matrix<Value>& occupation, io::Matrix<Value>& correlation)
        {
            jsx::value jOccupation;
            
            
            std::cout << "Calculating occupation ... " << std::flush;
            
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
        
            
            jsx::value jn = jsx::array_t(jHybMatrix.size());
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                linalg::mult<Value>('c', 'n', 1., jParams("operators")(i), jParams("operators")(i), .0, jn(i));
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    std::string entry = jHybMatrix(i)(j).string();
                    
                    if(!entry.empty()) {
                        jsx::value jOccupation;
                        linalg::mult<Value>('c', 'n', 1., jParams("operators")(i), jParams("operators")(j), .0, jOccupation);
                        occupation(i, j) = linalg::trace<Value>(jDensityMatrix, jOccupation);
                    }
                    
                    jsx::value jCorrelation;
                    linalg::mult<Value>('n', 'n', 1., jn(i), jn(j), .0, jCorrelation);
                    correlation(i, j) = linalg::trace<Value>(jDensityMatrix, jCorrelation);
                }
            
            auto occupation_entries = func::get_entries(occupation, jHybMatrix);
            
            for(auto const& entry : occupation_entries)
                jOccupation[entry.first] = io::Vector<Value>{{entry.second}};
            
            occupation = func::get_matrix(occupation_entries, jHybMatrix);
            
            std::cout << "Ok" << std::endl;
            
            
            return jOccupation;
        }
        
    }
}


#endif









