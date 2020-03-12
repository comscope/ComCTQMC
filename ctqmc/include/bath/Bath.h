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
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = default;
        Matrix& operator=(Matrix const&) = delete;
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

    
    template<typename Value>
    struct Bath {
		Bath() : det_(1.), B_(0) {};
        Bath(Bath const& bath) = delete;
        Bath(Bath&& bath) = delete;
        Bath& operator=(Bath const& bath) = delete;
        Bath& operator=(Bath&& bath) = delete;
        ~Bath() = default;

		std::vector<Operator<Value>> const& opsL() const { return opsL_;};
		std::vector<Operator<Value>> const& opsR() const { return opsR_;};
        
        void insertL(ut::KeyType key, int flavor) {
            posL_[key] = opsL_.size(); opsL_.push_back(Operator<Value>(key, flavor));
        };
        void insertR(ut::KeyType key, int flavor) {
            posR_[key] = opsR_.size(); opsR_.push_back(Operator<Value>(key, flavor));
        };
		
        Matrix<Value>& B() { return B_;};
        Matrix<Value> const& B() const { return B_;};
        
        Value sign() const { return det_/std::abs(det_);};
	
        void clean(Hyb<Value> const& hyb) {
            if(opsL_.size() != opsR_.size())
                throw std::runtime_error("bath::bath::clean: invalid configuration.");
            
            int const N = opsL_.size();
            
            det_ = 1.;  if(N != B_.dim()) B_ = Matrix<Value>(N);  

			if(N) {
				Matrix<Value> toInvert(N);

				for(int j = 0; j < N; ++j) 					
					for(int i = 0; i < N; ++i) 
						toInvert.at(i, j) = hyb(opsL_[i].flavor(), opsR_[j].flavor(), opsL_[i].key() - opsR_[j].key());
				
                std::fill(B_.data(), B_.data() + N*N, .0);
                for(int i = 0; i < N; ++i) B_.at(i, i) = 1.;
				
				int* ipiv = new int[N]; int info;
				gesv(&N, &N, toInvert.data(), &N, ipiv, B_.data(), &N, &info);
				for(int i = 0; i < N; ++i) 
					det_ *= (ipiv[i] != i + 1 ? -toInvert.at(i, i) : toInvert.at(i, i));
				delete[] ipiv;
			}
		}
        
        template<typename Type>
        void add(Hyb<Value> const& hyb, Type type) {
            if(update_.get() == nullptr)
                update_.reset(new Update<Value, Type>());
            else if(update_->type().val_ != get_type_id<Type>().val_)
                throw std::runtime_error("bath::add: update type missmatch");
            
            static_cast<Update<Value, Type>&>(*update_).add(type, *this, hyb);
        };
        
        double ratio(Hyb<Value> const& hyb) {
            return update_.get() == nullptr ? 1. : update_->ratio(*this, hyb);
        };
        
        int accept(Hyb<Value> const& hyb) {
            int sign = 1;
            if(update_.get() != nullptr) {
                sign = update_->accept(*this, hyb);  update_.reset(nullptr);
            }
            return sign;
        };
        
        void reject(Hyb<Value> const& hyb) {
            if(update_.get() != nullptr) {
                update_->reject(*this, hyb);  update_.reset(nullptr);
            }
        };

	private:
		Value det_;
		Matrix<Value> B_;
        
        std::vector<Operator<Value>> opsL_;
        std::vector<Operator<Value>> opsR_;
        std::map<ut::KeyType, int> posL_;
        std::map<ut::KeyType, int> posR_;
        
        std::unique_ptr<itf::Update<Value>> update_;
        
        template<typename, typename> friend struct Update;
	};
}

#endif
