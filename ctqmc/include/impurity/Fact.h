#ifndef CTQMC_INCLUDE_IMPURITY_FACT_H
#define CTQMC_INCLUDE_IMPURITY_FACT_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../Utilities.h"


namespace imp {

    template<typename Value>
    struct Fact {
        Fact() = default;
        Fact(Fact const&) = delete;
        Fact(Fact&&) = delete;
        Fact& operator=(Fact const&) = delete;
        ~Fact() = default;

        void operator*=(Value arg);
        
        void operator/=(Value arg);
        
        double ratio() const;
        
        void accept();
        
        void reject();
        
        Value sign() const;
        
    private:
        Value fact_{1.};
        Value factTry_{1.};
        
    };
    
}

#include "Fact.impl.h"

#endif
