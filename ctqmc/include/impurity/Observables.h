#ifndef CTQMC_IMPURITY_OBSERVABLES_H
#define CTQMC_IMPURITY_OBSERVABLES_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Algebra.h"
#include "Diagonal.h"
#include "Operators.h"
#include "../../../include/JsonX.h"

namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct BullaOperators {
            virtual int flavors() const = 0;
            virtual ~BullaOperators() = default;
        };
        
        template<typename Value>
        struct Occupation {
            virtual int flavors() const = 0;
            virtual ~Occupation() = default;
        };
        
        template<typename Value>
        struct BullaOccupation {
            virtual int flavors() const = 0;
            virtual ~BullaOccupation() = default;
        };
        
    };
    
    
    template<typename Mode, typename Value>
    struct BullaOperators : itf::BullaOperators<Value> {
        BullaOperators() = delete;
        BullaOperators(jsx::value const& jMPI, jsx::value const& jInteraction, jsx::value const& jOperators, itf::EigenValues const& eig);
        BullaOperators(BullaOperators const&) = delete;
        BullaOperators(BullaOperators&&) = delete;
        BullaOperators& operator=(BullaOperators const&) = delete;
        BullaOperators& operator=(BullaOperators&&) = delete;
        ~BullaOperators();
        
        int flavors() const {
            return flavors_;
        };
        Operator<Mode, Value> const& at(int f) const {
            return ops_[f];
        };
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> BullaOperators<Mode, Value>& get(itf::BullaOperators<Value>& bullaOpsItf) {
        return static_cast<BullaOperators<Mode, Value>&>(bullaOpsItf);
    };
    
    template<typename Mode, typename Value> BullaOperators<Mode, Value> const& get(itf::BullaOperators<Value> const& bullaOpsItf) {
        return static_cast<BullaOperators<Mode, Value> const&>(bullaOpsItf);
    };
    
    
    template<typename Mode, typename Value>
    struct Occupation : itf::Occupation<Value> {
        Occupation() = delete;
        Occupation(jsx::value const& jOperators, itf::EigenValues const& eig);
        Occupation(Occupation const&) = delete;
        Occupation(Occupation&&) = delete;
        Occupation& operator=(Occupation const&) = delete;
        Occupation& operator=(Occupation&&) = delete;
        ~Occupation();
        
        int flavors() const { return flavors_;};
        Operator<Mode, Value> const& at(int f) const { return ops_[f];};
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> Occupation<Mode, Value>& get(itf::Occupation<Value>& occItf) {
        return static_cast<Occupation<Mode, Value>&>(occItf);
    };
    
    template<typename Mode, typename Value> Occupation<Mode, Value> const& get(itf::Occupation<Value> const& occItf) {
        return static_cast<Occupation<Mode, Value> const&>(occItf);
    };
    
    
    template<typename Mode, typename Value>
    struct BullaOccupation : itf::BullaOccupation<Value> {
        BullaOccupation() = delete;
        BullaOccupation(jsx::value const& jParams, std::vector<double> const& filling, jsx::value jEigenValues, jsx::value const& jOperators, itf::EigenValues const& eig);
        BullaOccupation(BullaOccupation const&) = delete;
        BullaOccupation(BullaOccupation&&) = delete;
        BullaOccupation& operator=(BullaOccupation const&) = delete;
        BullaOccupation& operator=(BullaOccupation&&) = delete;
        ~BullaOccupation();
        
        int flavors() const { return flavors_;};
        Operator<Mode, Value> const& at(int f) const { return ops_[f];};
        
    private:
        int const flavors_;
        Operator<Mode, Value>* ops_;
    };
    
    template<typename Mode, typename Value> BullaOccupation<Mode, Value>& get(itf::BullaOccupation<Value>& bullaOccItf) {
        return static_cast<BullaOccupation<Mode, Value>&>(bullaOccItf);
    };
    
    template<typename Mode, typename Value> BullaOccupation<Mode, Value> const& get(itf::BullaOccupation<Value> const& bullaOccItf) {
        return static_cast<BullaOccupation<Mode, Value> const&>(bullaOccItf);
    };
    
}

#include "Observables.h"

#endif
