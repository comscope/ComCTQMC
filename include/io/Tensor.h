#ifndef INCLUDE_IO_TENSOR_H
#define INCLUDE_IO_TENSOR_H

#include <vector>
#include <complex>

#include "Vector.h"
#include "../JsonX.h"
#include "../../ctqmc/include/Utilities.h"

//scheiss b64 member variable !! Kack loesig im moment ...

namespace io {
    
    template<typename T> struct Tensor {
        inline static std::string name() { return name(T());};
        
        Tensor() = default;
        Tensor(std::size_t I, std::size_t J, std::size_t K, std::size_t L, T value = .0) : I_(I), J_(J), K_(K), L_(L), data_(I_*J_*K_*L_, value) {};
        Tensor(Tensor const&) = default;
        Tensor(Tensor&& other) noexcept : I_(other.I_), J_(other.J_),K_(other.K_), L_(other.L_), data_(std::move(other.data_)) { other.I_ = other.J_ = other.K_ = other.L_ = 0;};
        Tensor& operator=(Tensor const&) = default;
        Tensor& operator=(Tensor&& other) { I_ = other.I_; J_ = other.J_; K_ = other.K_; L_ = other.L_; other.I_ = other.J_ = other.K_ = other.L_ = 0; data_ = std::move(other.data_);  return *this;};
        ~Tensor() = default;

        int const& I() const { return I_;};
        int const& J() const { return J_;};
        int const& K() const { return K_;};
        int const& L() const { return L_;};
        
        T* data() { return data_.data();};
        T const* data() const { return data_.data();};
        
        T& operator()(int i, int j, int k, int l) { return data_.at(i + j*I_ + k*I_*J_ + l*I_*J_*K_);};
        T const& operator()(int i, int j, int k, int l) const { return data_.at(i + j*I_ + k*I_*J_ + l*I_*J_*K_);};
        
        Tensor& resize(int I, int J, int K, int L, T value = .0) {
            I_ = I; J_ = J; K_ = K; L_ = L; data_.resize(I_*J_*K_*L_, value);
            return *this;
        };
        Tensor& conj() {
            Tensor temp; temp.data_.resize(I_*J_*K_*L_);
            temp.I_ = J_; temp.J_ = I_; temp.K_ = K_; temp.L_ = L_;
            
            for(int i = 0; i < I_; ++i)
                for(int j = 0; j < J_; ++j)
                    for(int k = 0; k < K_; ++k)
                        for(int l = 0; l < L_; ++l)
                            temp(l, k, j, i) = ut::conj(operator()(i, j, k, l));
            
            return *this = std::move(temp);
        };
        
        void read(jsx::value const& source) {
            I_ = source(0).int64();
            J_ = source(1).int64();
            K_ = source(2).int64();
            L_ = source(3).int64();
            data_.read(source(4));
        };
        void write(jsx::value& dest) const {
            dest = jsx::array_t(5);
            
            dest(0) = static_cast<jsx::int64_t>(I());
            dest(1) = static_cast<jsx::int64_t>(J());
            dest(2) = static_cast<jsx::int64_t>(K());
            dest(3) = static_cast<jsx::int64_t>(L());
            data_.write(dest(4));
        };
        bool& b64() const {
            return data_.b64();
        };
        
    private:
        int I_ = 0, J_ = 0, K_ = 0, L_ = 0;;
        Vector<T> data_;

        inline static std::string name(double const&) { return "io::rtens";};
        inline static std::string name(std::complex<double> const&) { return "io::ctens";};
    };

    
    typedef Tensor<double> rtens;
    typedef Tensor<std::complex<double>> ctens;
    
    
};

#endif 
