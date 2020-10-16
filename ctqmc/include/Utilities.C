
#include "Utilities.h"

namespace ut {
    
    Zahl<double>::Zahl() : mantissa_(.0), exponent_(std::numeric_limits<int>::min()) {};
    
    Zahl<double>::Zahl(double x, double y) {   // constructs x*exp(y)
        if(std::isfinite(x) && std::isfinite(y)) {
            mantissa_ = std::frexp(x*std::exp(y - M_LN2*(static_cast<int>(y/M_LN2) + 1)), &exponent_);
            mantissa_ != .0 ? exponent_ += static_cast<int>(y/M_LN2) + 1 : exponent_ = std::numeric_limits<int>::min();
        } else
            throw std::runtime_error("ut::Zahl<double>: argument of constructor is not a number");
    };
    
    Zahl<double>& Zahl<double>::operator+=(Zahl const& arg) {
        int exp;
        if(exponent_ > arg.exponent_) {
            mantissa_ = std::frexp(mantissa_ + std::ldexp(arg.mantissa_, arg.exponent_ - exponent_), &exp);
            mantissa_ != .0 ? exponent_ += exp : exponent_ = std::numeric_limits<int>::min();
        } else {
            mantissa_ = std::frexp(arg.mantissa_ + std::ldexp(mantissa_, exponent_ - arg.exponent_), &exp);
            mantissa_ != .0 ? exponent_ = arg.exponent_ + exp : exponent_ = std::numeric_limits<int>::min();
        }
        return *this;
    };
    
    Zahl<double>& Zahl<double>::operator*=(Zahl const& arg) {
        int exp; mantissa_ = std::frexp(mantissa_*arg.mantissa_, &exp);
        mantissa_ != .0 ? exponent_ += (arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
        return *this;
    };
    
    Zahl<double>& Zahl<double>::operator/=(Zahl const& arg) {
        if(arg.mantissa_ == .0) throw std::runtime_error("ut::Zahl<double>: division by zero");
        
        int exp; mantissa_ = std::frexp(mantissa_/arg.mantissa_, &exp);
        mantissa_ != .0 ? exponent_ += (-arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
        return *this;
    };
    
    double Zahl<double>::get() const {
        return std::ldexp(mantissa_, exponent_);
    };
    
    
    Zahl<double> Zahl<double>::abs() const {
        Zahl temp = *this; temp.mantissa_ = std::abs(temp.mantissa_); return temp;
    };
    
    
    Zahl<complex>::Zahl() : mantissa_(.0), exponent_(std::numeric_limits<int>::min()) {};
    
    
    Zahl<complex>::Zahl(complex z, double y) {
        if(std::isfinite(z.real()) && std::isfinite(z.imag()) && std::isfinite(y)) {  //better this way than using constructor of Zahl<double> because of throw
            double const x = std::abs(z);
            mantissa_ = x != .0 ? std::frexp(x*std::exp(y - M_LN2*(static_cast<int>(y/M_LN2) + 1)), &exponent_)*(z/x) : .0;
            mantissa_ != .0 ? exponent_ += static_cast<int>(y/M_LN2) + 1 : exponent_ = std::numeric_limits<int>::min();
        } else
            throw std::runtime_error("ut::Zahl<complex>: argument of constructor is not a number");
    };
    
    
    Zahl<complex>& Zahl<complex>::operator+=(Zahl const& arg) {
        int exp;
        if(exponent_ > arg.exponent_) {
            auto const z = mantissa_ + (arg.mantissa_ != .0 ? arg.mantissa_*std::ldexp(1., arg.exponent_ - exponent_) : .0);
            auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
            mantissa_ != .0 ? exponent_ += exp : exponent_ = std::numeric_limits<int>::min();
        } else {
            auto const z = arg.mantissa_ + (mantissa_ != .0 ? mantissa_*std::ldexp(1., exponent_ - arg.exponent_) : .0);
            auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
            mantissa_ != .0 ? exponent_ = arg.exponent_ + exp : exponent_ = std::numeric_limits<int>::min();
        }
        return *this;
    };
    
    
    Zahl<complex>& Zahl<complex>::operator*=(Zahl const& arg) {
        auto const z = mantissa_*arg.mantissa_; int exp;
        auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
        mantissa_ != .0 ? exponent_ += (arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
        return *this;
    };
    
    
    Zahl<complex>& Zahl<complex>::operator/=(Zahl const& arg) {
        if(arg.mantissa_ == .0) throw std::runtime_error("ut::Zahl<complex>: division by zero");
        
        auto const z = mantissa_/arg.mantissa_; int exp;
        auto const x = std::abs(z); mantissa_ = x != .0 ? std::frexp(x, &exp)*(z/x) : .0;
        mantissa_ != .0 ? exponent_ += (-arg.exponent_ + exp) : exponent_ = std::numeric_limits<int>::min();
        return *this;
    };
    
    
    complex Zahl<complex>::get() const {
        return mantissa_*std::ldexp(1., exponent_);
    };
    
    
    Zahl<double> Zahl<complex>::abs() const {
        Zahl<double> temp; temp.mantissa_ = std::abs(mantissa_); temp.exponent_ = exponent_; return temp;
    };
    
}

