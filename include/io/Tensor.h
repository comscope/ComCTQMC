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
        Tensor(std::size_t I, std::size_t J, std::size_t K, std::size_t L);
        Tensor(Tensor const&) = default;
        Tensor(Tensor&& other) noexcept;
        Tensor& operator=(Tensor const&) = default;
        Tensor& operator=(Tensor&& other);
        ~Tensor() = default;

        int const& I() const { return I_;};
        int const& J() const { return J_;};
        int const& K() const { return K_;};
        int const& L() const { return L_;};
        
        T* data() { return data_.data();};
        T const* data() const { return data_.data();};
        std::vector<std::vector<int>> const& ijkl() const { return ijkl_;}
        
        //Will complatin if an element is not there
        T& operator()(int const i, int const j, int const k, int const l){ return this->operator()(this->index(i,j,k,l)); }
        T const& operator()(int const i, int const j, int const k, int const l) const { return this->operator()(this->index(i,j,k,l)); }
        
        //Use at to retrieve elements that may or may not be there (= 0 if not there)
        T const at(int const i, int const j, int const k, int const l) const { return this->at(this->index(i,j,k,l)); }
        
        std::string const& entry(int const i, int const j, int const k, int const l) const { return this->entry(this->index(i,j,k,l)); }
        bool const is(int const i, int const j, int const k, int const l) const { return this->is(this->index(i,j,k,l)); }
        
        void clear(){ this->data_.clear(); this->entries_.clear(); I_=0; J_=0; K_=0; L_=0;}
        
        
        inline T& operator()(int const i){
            if (data_.find(i) == data_.end())
                throw std::runtime_error("Cannot update element of " + this->name() + " that does not yet exist\n");
            return data_.at(i);
        };
        inline T const& operator()(int const i) const {
            if (data_.find(i) != data_.end())
                return data_.at(i);
            else
                return zero;
        };
        
        inline T const at(int const i) const {
            if (data_.find(i) != data_.end())
                return data_.at(i);
            else
                return zero;
        };
        
        inline std::string const& entry(int const i) const {
            if (entries_.find(i) != entries_.end())
                return entries_.at(i);
            else
                return empty_entry;
        };
                              

        inline bool is(int const i) const { return (data_.find(i) == data_.end()) ? false : true; }

        void emplace (int const i, int const j, int const k, int const l, std::string const& entry, T const v);
        
        Tensor& resize(int const I, int const J, int const K, int const L, T value = .0);
        Tensor& conj();
        
    private:
        int I_ = 0, J_ = 0, K_ = 0, L_ = 0;
        std::vector<std::vector<int>> ijkl_;
        std::map<int,T> data_;
        std::map<int,std::string> entries_;

        T const zero = 0.;
        std::string const empty_entry = "";
        
        inline static std::string name(double const&) { return "io::rtens";};
        inline static std::string name(std::complex<double> const&) { return "io::ctens";};
        
        inline int index(int const i, int const j, int const k, int const l) const { return i + j*I_ + k*I_*J_ + l*I_*J_*K_;}
        inline void emplace (int const i, std::string const& entry, T const v){
            entries_.emplace(i,entry);
            data_.emplace(i,v);
        }
    };
    
    
    typedef Tensor<double> rtens;
    typedef Tensor<std::complex<double>> ctens;
    
    
};

#include "Tensor.impl.h"

#endif 
