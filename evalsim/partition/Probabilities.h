#ifndef EVALSIM_PARTITION_PROBABILITIES_H
#define EVALSIM_PARTITION_PROBABILITIES_H


#include "ReadDensityMatrix.h"
#include "ReadHamiltonian.h"


#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/options/Options.h"
#include "../../include/atomic/Generate.h"

namespace evalsim {
    
    namespace partition {
        
        //Less than for vectors
        struct VecLess {
            VecLess(std::size_t size, std::size_t start, std::size_t end);
            bool operator()(std::vector<double> const& lhs, std::vector<double> const& rhs) const;
        private:
            std::size_t const size_, start_, end_;
        };
        
        
        double truncate(double val, int prec);
        
        ut::complex truncate(ut::complex val, int prec);

        template<typename Value>
        std::vector<std::string> get_surviving_qn(jsx::value const& jPartition);
        
        template<typename Value>
        jsx::array_t get_eigenstate_probabilities(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
    
        template<typename Value>
        jsx::array_t get_eigenstate_energies(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
    
        template<typename Value>
        jsx::array_t get_eigenstate_qn(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
        
        
        template<typename Value>
        jsx::value get_eigenstates(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
        
        double serial_number(int sector, int i);
        
        template<typename Value>
        jsx::value get_probabilities(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
   
    }
    
}


#endif









