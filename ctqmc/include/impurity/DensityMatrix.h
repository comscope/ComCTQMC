#ifndef CTQMC_INCLUDE_IMPURITY_DENSITYMATRIX_H
#define CTQMC_INCLUDE_IMPURITY_DENSITYMATRIX_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>

#include "Algebra.h"
#include "Diagonal.h"
#include "Operators.h"
#include "Product.h"
#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct DensityMatrix {
            virtual ut::Flag surviving(itf::EigenValues const&) = 0;
            virtual ut::Flag decide(ut::Zahl<double> const&, itf::Product<Value>&, itf::Batcher<Value>&) = 0;
            virtual std::vector<int>::const_iterator begin() const = 0;
            virtual std::vector<int>::const_iterator end() const = 0;
            virtual ut::Zahl<Value> Z() const = 0;
            virtual Value weight(int) const = 0;
            virtual Value sign() const = 0;
            virtual ~DensityMatrix() = default;
        };
        
    };
    
    struct Bound {
        int sec; double ln; ut::Zahl<double> value;
    };
    
    inline int compare(Bound const& a, Bound const& b) {
        return a.ln > b.ln;
    };
    
    //-----------------------------------------------------------------DENSITYMATRIX---------------------------------------------------------------------------------
    template<typename Mode, typename Value>
    struct DensityMatrix : itf::DensityMatrix<Value> {
        typedef std::vector<int>::const_iterator iterator;
        
        DensityMatrix() = default;
        DensityMatrix(itf::Product<Value>& product, itf::EigenValues const& eig);
        DensityMatrix(DensityMatrix const&) = delete;
        DensityMatrix(DensityMatrix&&) = default;
        DensityMatrix& operator=(DensityMatrix const&) = delete;
        DensityMatrix& operator=(DensityMatrix&&) = default;
        ~DensityMatrix() = default;
        
        ut::Flag surviving(itf::EigenValues const& eig);

        ut::Flag decide(ut::Zahl<double> const& thresh, itf::Product<Value>& product, itf::Batcher<Value>& batcher);
        
        std::vector<int>::const_iterator begin() const { return sectors_.begin();};
        std::vector<int>::const_iterator end() const { return sectors_.end();};
        
        Matrix<Mode, Value> const& mat(int s) const { return static_cast<Operator<Mode, Value> const*>(op_.get())->mat(s);};
        
        ut::Zahl<Value> Z() const { return Z_;};
        Value weight(int s) const { return (z_[s]/Z_).get();};
        
        Value sign() const;
        
    private:
        std::unique_ptr<Operator<Mode, Value>> op_;
        
        int level_;
        std::vector<int> sectors_;
        std::vector<Bound> bounds_;
        std::vector<Bound>::iterator bound_;
        
        ut::Zahl<Value> Z_; std::vector<ut::Zahl<Value>> z_;
    };
    
    template<typename Mode, typename Value> DensityMatrix<Mode, Value>& get(itf::DensityMatrix<Value>& densityMatrixItf) {
        return static_cast<DensityMatrix<Mode, Value>&>(densityMatrixItf);
    };
    
    template<typename Mode, typename Value> DensityMatrix<Mode, Value> const& get(itf::DensityMatrix<Value> const& densityMatrixItf) {
        return static_cast<DensityMatrix<Mode, Value> const&>(densityMatrixItf);
    };
    
}

#include "DensityMatrix.impl.h"

#endif
