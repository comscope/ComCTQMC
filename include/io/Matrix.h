#ifndef INCLUDE_IO_MATRIX_H
#define INCLUDE_IO_MATRIX_H

#include <vector>
#include <complex>

#include "Vector.h"
#include "../JsonX.h"
#include "../../ctqmc/include/Utilities.h"

//scheiss b64 member variable !! Kack loesig im moment ...

namespace io {
    
    template<typename T> struct Matrix {
        inline static std::string name() { return name(T());};
        
        Matrix() = default;
        Matrix(std::size_t I, std::size_t J, T value = .0);
        Matrix(Matrix const&) = default;
        Matrix(Matrix&& other) noexcept;
        Matrix& operator=(Matrix const&) = default;
        Matrix& operator=(Matrix&& other);
        ~Matrix() = default;

        int const& I() const { return I_;};
        int const& J() const { return J_;};
        
        T* data() { return data_.data();};
        T const* data() const { return data_.data();};
        
        T& operator()(int i, int j) { return data_.at(i + j*I_);};
        T const& operator()(int i, int j) const { return data_.at(i + j*I_);};
        
        Matrix& resize(int I, int J, T value = .0);
        Matrix& conj();
        
        Matrix conj() const;
        
        void read(jsx::value const& source);
        void write(jsx::value& dest) const;
        
        bool& b64() const {
            return data_.b64();
        };
        
    private:
        int I_ = 0, J_ = 0;
        Vector<T> data_;

        inline static std::string name(double const&) { return "io::rmat";};
        inline static std::string name(std::complex<double> const&) { return "io::cmat";};
    };

    
    typedef Matrix<double> rmat;
    typedef Matrix<std::complex<double>> cmat;
    
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Patch    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    template<typename T> struct PrettyMatrix {
        inline static std::string name() { return name(T());};
        
        PrettyMatrix() = default;
        PrettyMatrix(std::size_t I, std::size_t J, T value = .0);
        PrettyMatrix(PrettyMatrix const&) = default;
        PrettyMatrix(PrettyMatrix&& other) noexcept;
        PrettyMatrix& operator=(PrettyMatrix const&) = default;
        PrettyMatrix& operator=(PrettyMatrix&& other);
        ~PrettyMatrix() = default;
        
        int const& I() const { return I_;};
        int const& J() const { return J_;};
        
        T& operator()(int i, int j) { return data_.at(J_*i + j);};
        T const& operator()(int i, int j) const { return data_.at(J_*i + j);};
        
        PrettyMatrix& resize(int I, int J, T value = .0);

        void read(jsx::value const& source);
        void write(jsx::value& dest) const;
        
    private:
        int I_ = 0, J_ = 0;
        std::vector<T> data_;
        
        inline static std::string name(double const&) { return "io::prettyrmat";};
        inline static std::string name(std::complex<double> const&) { return "io::prettycmat";};
    };
    
    
    
    typedef PrettyMatrix<double> prettyrmat;
    typedef PrettyMatrix<std::complex<double>> prettycmat;
};

#include "Matrix.impl.h"

#endif 
