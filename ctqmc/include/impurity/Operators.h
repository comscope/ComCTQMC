#ifndef CTQMC_INCLUDE_IMPURITY_OPERATORS_H
#define CTQMC_INCLUDE_IMPURITY_OPERATORS_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>


#include "Algebra.h"
#include "BitSet.h"
#include "Diagonal.h"
#include "Utilities.h"
#include "../Utilities.h"
#include "../../../include/linalg/LinAlg.h"
#include "../../../include/linalg/Operators.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct Operator {
            virtual SectorNorm const& map(int) const = 0;
            virtual ~Operator() = default;
        };
        
        template<typename Value>
        struct Operators {
            virtual int flavors() const = 0;
            virtual Operator<Value> const& at(int) const = 0;
            virtual ~Operators() = default;
        };
        
    };
    
    //-----------------------------------------------------------------OPERATOR--------------------------------------------------------------------------------
    template<typename Mode, typename Value>
    struct Operator : itf::Operator<Value> {
        Operator() = delete;
        Operator(itf::EigenValues const& eigItf);
        
        Operator(char const option, itf::EigenValues const& eigItf);
        
        //Supports parallelization by pre-computing norms
        //otherwise, if no norms are passed, each rank computes all norms
        Operator(jsx::value const& jOperator, itf::EigenValues const& eigItf, io::rvec const& norms = {});
        Operator(Operator const&) = delete;
        Operator(Operator&&) = delete;
        Operator& operator=(Operator const&) = delete;
        Operator& operator=(Operator&&) = delete;
        ~Operator();
        
        int isMap(int s) const;
        SectorNorm& set_map(int s);
        
        SectorNorm& map(int s);
        SectorNorm const& map(int s) const;
        
        int isMat(int s) const;
        template<typename... Args> Matrix<Mode, Value>& mat(int s, Args&& ... args);
        
        Matrix<Mode, Value>& mat(int s);
        Matrix<Mode, Value> const& mat(int s) const;
        
        void assign(std::vector<int> const& sectors, SectorNormPtrs& arg);
        int missing(SectorNormPtrs& missing, SectorNormPtrs const& requested);
        void map(SectorNormPtrs& arg) const;
        
    private:
        EigenValues<Mode> const& eig_;
        
        BitSet isMap_, isMat_;
        SectorNorm* const map_;
        Matrix<Mode, Value>* const mat_;
    };
    
    template<typename Mode, typename Value> Operator<Mode, Value>& get(itf::Operator<Value>& operatorItf) {
        return static_cast<Operator<Mode, Value>&>(operatorItf);
    };
    
    template<typename Mode, typename Value> Operator<Mode, Value> const& get(itf::Operator<Value> const& operatorItf) {
        return static_cast<Operator<Mode, Value> const&>(operatorItf);
    };
    
    
    
    template<typename Mode, typename Value>
    struct Operators : itf::Operators<Value> {
        Operators() = delete;
        Operators(jsx::value const& jParams, jsx::value const& jOperators, itf::EigenValues const& eigItf);
        
        Operators(Operators const&) = delete;
        Operators(Operators&&) = delete;
        Operators& operator=(Operators const&) = delete;
        Operators& operator=(Operators&&) = delete;
        ~Operators();
        
        int flavors() const { return flavors_;};
        
        itf::Operator<Value> const& at(int i) const { return ops_[i];};
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> Operators<Mode, Value>& get(itf::Operators<Value>& operatorsItf) {
        return static_cast<Operators<Mode, Value>&>(operatorsItf);
    };
    
    template<typename Mode, typename Value> Operators<Mode, Value> const& get(itf::Operators<Value> const& operatorsItf) {
        return static_cast<Operators<Mode, Value> const&>(operatorsItf);
    };
    
}

#include "Operators.impl.h"

#endif  
