#ifndef CTQMC_INCLUDE_BATH_HYB_H
#define CTQMC_INCLUDE_BATH_HYB_H

#include <cmath>
#include <iostream>
#include <sstream>

#include "../Utilities.h"
#include "../../../include/atomic/Generate.h"

namespace bath {
    
    struct Block {
        std::vector<int>& flavorsL() { return flavorsL_;};
        std::vector<int>& flavorsR() { return flavorsR_;};
        
        std::vector<int> const& flavorsL() const { return flavorsL_;};
        std::vector<int> const& flavorsR() const { return flavorsR_;};
        
    private:
        std::vector<int> flavorsL_;
        std::vector<int> flavorsR_;
    };
    

    template<typename Value>
    struct Fit {
        Fit(double const beta, std::vector<std::complex<double>> const& hyb, std::vector<std::complex<double>> const& hybTransp);
        
        std::size_t N() const { return N_;};
        Value moment() const { return moment_;};
        Value eps() const { return eps_;};
    private:
        std::size_t N_;
        Value moment_;
        Value eps_;
        
        void check(ut::complex, ut::complex, ut::complex) {};
        void check(ut::complex pos, ut::complex neg, double);
    };
    
    
    template<typename Value>
    struct Simple {
        Simple() = delete;
        Simple(jsx::value const& jParams, std::vector<std::complex<double>> const& hyb, std::vector<std::complex<double>> hybTransp);
        Simple(Simple const&) = delete;
        Simple(Simple&&) = default;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
        
        Value get(ut::KeyType key) const;
        
        ~Simple() = default;
    private:
        std::size_t const I_;
        double const fact_;
        
        std::vector<Value> data_;
    };
    
    
    template<typename Value>
    struct Hyb {
        Hyb() = delete;
        Hyb(jsx::value const& jParams, jsx::value const& jMatrix, jsx::value jFunctions);
        Hyb(Hyb const&) = delete;
        Hyb(Hyb&&) = delete;
        Hyb& operator=(Hyb const& hyb) = delete;
        Hyb& operator=(Hyb&&) = delete;
        ~Hyb() = default;
        
        std::vector<Block> const& blocks() const { return blocks_;};
        std::size_t bath(int flavor) const { return bath_[flavor];};
        bool isR(int flavor) const { return isR_[flavor];};
        
        Value operator()(int flavorL, int flavorR, ut::KeyType key) const;
        
        int flavors() const { return flavors_; }
        
    private:
        int const flavors_;
        std::map<std::string, Simple<Value>> data_;
        std::vector<Simple<Value> const*> matrix_;

        std::vector<Block> blocks_;
        std::vector<int> bath_;
        std::vector<bool> isR_;
    };
}

#include "Hyb.impl.h"

#endif
