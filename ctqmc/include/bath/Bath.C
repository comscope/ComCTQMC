
#include "Bath.h"

//Vielleicht insert und erase noch verändern wie in CTBI hürütütütüt pluschter pluschter

namespace bath {
    
    template<typename Value>
    Bath<Value>::Bath(Bath const& bath) :
    opsL_(bath.opsL_), opsR_(bath.opsR_),
    posL_(bath.posL_), posR_(bath.posR_),
    det_(bath.det_), B_(bath.B_){};
    
    template<typename Value>
    Bath<Value>::Bath(Bath<Value>&& bath) :
    opsL_(std::move(bath.opsL_)), opsR_(std::move(bath.opsR_)),
    posL_(std::move(bath.posL_)), posR_(std::move(bath.posR_)),
    det_(std::move(bath.det_)), B_(std::move(bath.B_)){};
    
    template<typename Value>
    Bath<Value>& Bath<Value>::operator=(Bath const& bath){
        opsL_ = bath.opsL_; opsR_ = bath.opsR_;
        posL_ = bath.posL_; posR_ = bath.posR_;
        det_ = bath.det_; B_ = bath.B_;
        return *this;
        
    }
    
    template<typename Value>
    Bath<Value>& Bath<Value>::operator=(Bath&& bath){
        opsL_ = std::move(bath.opsL_); opsR_ = std::move(bath.opsR_);
        posL_ = std::move(bath.posL_); posR_ = std::move(bath.posR_);
        det_ = std::move(bath.det_); B_ = std::move(bath.B_);
        return *this;
    }
    
    
    template<typename Value>
    void Bath<Value>::insertL(ut::KeyType key, int flavor) {
        posL_[key] = opsL_.size(); opsL_.push_back(Operator<Value>(key, flavor));
    };
    
    template<typename Value>
    void Bath<Value>::insertR(ut::KeyType key, int flavor) {
        posR_[key] = opsR_.size(); opsR_.push_back(Operator<Value>(key, flavor));
    };
    
    template<typename Value>
    int Bath<Value>::eraseL(ut::KeyType key) {
        int sign = 1;
        auto const itL = posL_.find(key); auto const posL = itL->second;
        if(posL != opsL_.size() - 1)
            sign*=-1;
        posL_[opsL_.back().key()] = posL; opsL_[posL] = opsL_.back();
        posL_.erase(itL); opsL_.pop_back();
        return sign;
    };
    
    template<typename Value>
    int Bath<Value>::eraseR(ut::KeyType key) {
        int sign = 1;
        auto const itR = posR_.find(key); auto const posR = itR->second;
        if(posR != opsR_.size() - 1)
            sign*=-1;
        posR_[opsR_.back().key()] = posR; opsR_[posR] = opsR_.back();
        posR_.erase(itR); opsR_.pop_back();
        return sign;
    };
    
    
    template<typename Value>
    void Bath<Value>::clean(Hyb<Value> const& hyb) {
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
    
    template<typename Value>
    double Bath<Value>::ratio(Hyb<Value> const& hyb) {
        return update_.get() == nullptr ? 1. : update_->ratio(*this, hyb);
    };
    
    template<typename Value>
    int Bath<Value>::accept(Hyb<Value> const& hyb) {
        int sign = 1;
        if(update_.get() != nullptr) {
            sign = update_->accept(*this, hyb);  update_.reset(nullptr);
        }
        return sign;
    };
    
    template<typename Value>
    void Bath<Value>::reject(Hyb<Value> const& hyb) {
        if(update_.get() != nullptr) {
            update_->reject(*this, hyb);  update_.reset(nullptr);
        }
    };
		
    template struct Bath<double>;
    template struct Bath<ut::complex>;
        
}
