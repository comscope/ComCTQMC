#ifndef INCLUDE_LINALG_OPERATORS_H
#define INCLUDE_LINALG_OPERATORS_H

#include "../io/Vector.h"
#include "../io/Matrix.h"
#include "LinAlg.h"


namespace linalg {
    
    template<typename Value>
    jsx::value conj(jsx::value const& jOp);


    template<typename Value>
    jsx::value diag_to_operator(jsx::value const& jDiag);

    template<typename Value>
    void mult(char transA, char transB, double alpha, jsx::value const& jA, jsx::value const& jB, double beta, jsx::value& jC);


    template<typename Value>
    Value trace(jsx::value const& jA, jsx::value const& jB);

}

#include "Operators.impl.h"

#endif










