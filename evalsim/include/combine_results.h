#ifndef EVALSIM_INCLUDE_COMBINE_RESULTS
#define EVALSIM_INCLUDE_COMBINE_RESULTS

#include "../../include/io/Vector.h"
#include "../../include/JsonX.h"
#include "../../include/mpi/Utilities.h"

namespace evalsim{

struct MeasurementCombiner{
    
    MeasurementCombiner(){};
    
    void combine(jsx::value const& jParams, std::string const name);
    
    jsx::value& data() {return data_;}
    
    
    void add(jsx::value const& jIn, jsx::value & jMeas, double const factor, bool const first);
    
    
  
private:
    jsx::value data_;
};

}

#endif

