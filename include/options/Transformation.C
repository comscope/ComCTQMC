#include "Transformation.h"

namespace opt {
    
    
    template<typename Value>
    Transformation<Value>::Transformation(int const N, jsx::value jTransformation) :
    t_(jTransformation.is<jsx::empty_t>() ? io::PrettyMatrix<Value>(N, N) : jsx::at<io::PrettyMatrix<Value>>(jTransformation)) {
        if(jTransformation.is<jsx::empty_t>())
            for(int f = 0; f < N; ++f) t_(f, f) = 1.;
    }

    template<typename Value>
    Value Transformation<Value>::operator()(int i, int j) const {
        return t_(i, j);
    };
    template<typename Value>
    int Transformation<Value>::I() const { return t_.I();};

    template<typename Value>
    int Transformation<Value>::J() const { return t_.J();};

    
    template struct Transformation<double>;
    template struct Transformation<ut::complex>;
    
};



