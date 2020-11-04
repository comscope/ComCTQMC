#ifndef INCLUDE_OPTIONS_TRANSFORMATION_H
#define INCLUDE_OPTIONS_TRANSFORMATION_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <set>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>

#include "../JsonX.h"
#include "../io/Matrix.h"


namespace opt {
    
    template<typename Value>
    struct Transformation {
        Transformation(int const N, jsx::value jTransformation);
        
        Value operator()(int i, int j) const;
        
        int I() const;
        int J() const;
        
    private:
        io::PrettyMatrix<Value> t_;
    };
    
};

#include "Transformation.impl.h"
        
#endif






