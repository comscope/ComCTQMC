
#include "Fact.h"

namespace imp {

    template<typename Value>
    void Fact<Value>::operator*=(Value arg) {
        factTry_ *= arg;
    };
    
    template<typename Value>
    void Fact<Value>::operator/=(Value arg) {
        factTry_ /= arg;
    };
    
    template<typename Value>
    double Fact<Value>::ratio() const {
        return std::abs(factTry_/fact_);
    };
    
    template<typename Value>
    void Fact<Value>::accept() {
        fact_ = factTry_;
    };
    
    template<typename Value>
    void Fact<Value>::reject() {
        factTry_ = fact_;
    };
    
    template<typename Value>
    Value Fact<Value>::sign() const {
        return fact_/std::abs(fact_);
    };
    
    template struct Fact<double>;
    template struct Fact<ut::complex>;
    
}

