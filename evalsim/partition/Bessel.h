#ifndef EVALSIM_PARTITON_BESSEL_H
#define EVALSIM_PARTITON_BESSEL_H

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>

namespace evalsim {
    
    namespace be {
        
        struct Bessel {
            Bessel() = delete;
            Bessel(std::size_t const N, std::size_t const M);
            Bessel(Bessel const&) = delete;
            Bessel(Bessel&&) = default;
            Bessel& operator=(Bessel const&) = delete;
            Bessel& operator=(Bessel&&) = delete;
            ~Bessel() = default;
            
            double operator()(std::size_t n, std::size_t m) {
                return data_[n][m];
            };
            
        private:
            std::vector<std::vector<double>> data_;
            
            std::vector<double> bessel_forward(std::size_t const n_high, std::size_t const m);
            
            std::vector<double> bessel_backward(std::size_t const n_low, std::size_t const n_high, std::size_t const m, double const val_before_low);
            
            std::vector<double> bessel(std::size_t const n, std::size_t const m);
        };
        
    }
    
}

#endif


