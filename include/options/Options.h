#ifndef INCLUDE_OPTIONS_OPTIONS_H
#define INCLUDE_OPTIONS_OPTIONS_H

#include "Basis.h"
#include "Interaction.h"
#include "Observable.h"

#include "../JsonX.h"
#include "../io/Vector.h"
#include "../io/Matrix.h"
#include "../mpi/Utilities.h"
#include "../atomic/Generate.h"
#include "../../ctqmc/include/Utilities.h"

#include <ctime>


namespace opt {

    template<typename Value>
    jsx::value transform(Basis<Value> const& basis, Interaction const& interaction);

    template<typename Value>
    jsx::value truncate(io::Vector<Value> const& interaction, double const truncate);

    template<typename Value>
    void complete_hloc(jsx::value& jParams);

    template<typename Value>
    void complete_qn(jsx::value const& jParams, jsx::value& jqn);

    template<typename Value>
    jsx::value transform(Observable const& tensor, Transformation<Value> const& transformation, bool ising);

    template <typename Value>
    jsx::value read_observable(jsx::value const& jObservable, int const N, bool ising);

    template<typename Value>
    void complete_observables(jsx::value& jParams, jsx::value& jObservables, bool ising);
    
};

#include "Options.impl.h"

#endif //OPTIONS


