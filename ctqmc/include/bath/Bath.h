#ifndef CTQMC_INCLUDE_BATH_BATH_H
#define CTQMC_INCLUDE_BATH_BATH_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Hyb.h"
#include "../Utilities.h"
#include "../../../include/BlasLapack.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace bath {
    
    template<typename> struct Bath;
    

    struct TypeId {                                    // should be joined with JsonX.h ...
        using FuncPtr = TypeId(*)();
        TypeId(FuncPtr val) : val_(val) {};
    private:
        FuncPtr val_;
        template<typename> friend struct Bath;
    };
    
    template<typename T> TypeId get_type_id() { return &get_type_id<T>;};
    
    
    namespace itf {
        
        template<typename Value>
        struct Update {
            virtual TypeId type() = 0;
            virtual double ratio(Bath<Value> const&, Hyb<Value> const&) = 0;
            virtual int accept(Bath<Value>&, Hyb<Value> const&) = 0;
            virtual void reject(Bath<Value>&, Hyb<Value> const&) = 0;
            virtual ~Update() = default;
        };
        
    }
    
    template<typename, typename> struct Update;
    
    
    template<typename Value>
    struct Matrix {
        Matrix() = delete;
        explicit Matrix(int dim) : dim_(dim), data_(dim_*dim_) {};
        Matrix(Matrix const&) = default;
        Matrix(Matrix&&) = default;
        Matrix& operator=(Matrix const&) = default;
        Matrix& operator=(Matrix&&) = default;
        ~Matrix() = default;
        
        int dim() const { return dim_;};
        Value& at(int i, int j) { return data_[i + dim_*j];};
        Value const& at(int i, int j) const { return data_[i + dim_*j];};
        Value* data() { return data_.data();};
        Value const* data() const { return data_.data();};
        Value* data(int i, int j) { return data_.data() + i + j*dim_;};
        Value const* data(int i, int j) const { return data_.data() + i + j*dim_;};
        
    private:
        int dim_;
        std::vector<Value> data_;
    };
    
    
    template<typename Value>
    struct Operator {
        Operator(ut::KeyType key, int flavor) : key_(key), flavor_(flavor) {};
        
        ut::KeyType key() const { return key_;};
        int flavor() const { return flavor_;};
        Value& bulla() const { return bulla_;};
        
    private:
        ut::KeyType key_;
        int flavor_;
        mutable Value bulla_;
    };
/*
 std::vector<Operator<Value>> opsL_;
 std::vector<Operator<Value>> opsR_;
 std::map<ut::KeyType, int> posL_;
 std::map<ut::KeyType, int> posR_;
 */
    


    template<typename Value>
    struct Bath {
		Bath() : det_(1.), B_(0) {};
        
        //do I need to bring the update pointer along? seems unlikely.
        Bath(Bath const& bath);
        Bath(Bath&& bath);
        Bath& operator=(Bath const& bath);
        Bath& operator=(Bath&& bath);
        ~Bath() = default;

		std::vector<Operator<Value>> const& opsL() const { return opsL_;};
		std::vector<Operator<Value>> const& opsR() const { return opsR_;};
        
        void insertL(ut::KeyType key, int flavor);
        void insertR(ut::KeyType key, int flavor);
        
        int eraseL(ut::KeyType key);
        
        int eraseR(ut::KeyType key);
		
        Matrix<Value>& B() { return B_;};
        Matrix<Value> const& B() const { return B_;};
        
        Value sign() const { return det_/std::abs(det_);};
	
        void clean(Hyb<Value> const& hyb);
        
        template<typename Type>
        void add(Hyb<Value> const& hyb, Type type);
        
        double ratio(Hyb<Value> const& hyb);
        
        int accept(Hyb<Value> const& hyb);
        
        void reject(Hyb<Value> const& hyb);

	private:
        
        std::vector<Operator<Value>> opsL_;
        std::vector<Operator<Value>> opsR_;
        std::map<ut::KeyType, int> posL_;
        std::map<ut::KeyType, int> posR_;
        
        Value det_;
        Matrix<Value> B_;
        
        std::unique_ptr<itf::Update<Value>> update_;
        
        template<typename, typename> friend struct Update;
	};
}

#include "Bath.impl.h"

#endif
