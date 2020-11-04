#ifndef CTQMC_INCLUDE_IMPURITY_DYNAMIC_H
#define CTQMC_INCLUDE_IMPURITY_DYNAMIC_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


#include "../Utilities.h"
#include "../../../include/atomic/Generate.h"
#include "../../../include/mpi/Utilities.h"
#include "../../../include/JsonX.h"

//schä dömande a la lüünö si tü vulä ancor de mua e ell ma di wa tö fär encüle sche le gräk !! hahahhahahhahaha


// Todo: get quantum number change as a double which should be given by the operator

namespace imp {
    
    namespace itf {
        
        struct Dynamic {
            Dynamic() = default;
            Dynamic(Dynamic const&) = delete;
            Dynamic(Dynamic&&) = delete;
            Dynamic& operator=(Dynamic const&) = delete;
            Dynamic& operator=(Dynamic&&) = delete;
            virtual ~Dynamic() = default;
            
            virtual void insert(ut::KeyType key, int flavor) {};
            virtual void erase(ut::KeyType key, int flavor) {};
            
            virtual ut::Zahl<double> ratio() { return 1.;};
            
            virtual void accept() {};
            virtual void reject() {};
            
            virtual double energy() const { return .0;};
            virtual double qkinks(int qn) const { return .0;};
            virtual double fkinks(int flavor, ut::KeyType key) const { return .0;};
            
            virtual void clean() {};
        };
        
    };
    
    
    //TODO: check if interaction matrix is symmetric (because the interaction is assumed to be real in this implementation) !
    struct Simple {
        Simple() = delete;
        Simple(jsx::value const& jParams, jsx::value jDyn);
        Simple(Simple const&) = delete;
        Simple(Simple&&) = delete;
        Simple& operator=(Simple const&) = delete;
        Simple& operator=(Simple&&) = delete;
        ~Simple() = default;
        
        int size() const {
            return size_;
        };
        
        double shift(int sec) const {
            return shift_[sec];
        };
        
        double D0Qq(int sec, int qn) const {
            return D0Qq_[qn][sec];
        };
        
        double D0Qf(int sec, int flavor) const {
            return D0Qf_[flavor/2][sec];
        };
        
        double Lfq(int flavor, int qn, ut::KeyType key) const;
        
        double Lff(int flavorI, int flavorJ, ut::KeyType key) const;
        
        double Kff(int flavorI, int flavorJ, ut::KeyType key) const;

    private:
        int const size_;
        
        int nIt_;
        std::vector<double> shift_;
        std::vector<std::vector<double>> D0Qq_;
        std::vector<std::vector<double>> D0Qf_;
        std::vector<std::vector<std::vector<double>>> Lfq_;
        std::vector<std::vector<std::vector<double>>> Lff_;
        std::vector<std::vector<std::vector<double>>> Kff_;
    };

    
    struct Dynamic : itf::Dynamic {
        Dynamic() = delete;
        Dynamic(Simple const& func);
        Dynamic(Dynamic const&) = delete;
        Dynamic(Dynamic&&) = delete;
        Dynamic& operator=(Dynamic const&) = delete;
        Dynamic& operator=(Dynamic&&) = delete;
        ~Dynamic() = default;
        
        void insert(ut::KeyType key, int flavor);
        void erase(ut::KeyType key, int flavor);
        
        ut::Zahl<double> ratio();
        
        void accept();
        
        void reject();
        
        double energy() const;
        
        double qkinks(int qn) const;
        
        double fkinks(int flavor, ut::KeyType key) const;

        void clean();
    
    private:
        struct Entry {
            Entry() {};
            Entry(ut::KeyType key, int flavor, int state) : key(key), flavor(flavor), state(state) {};
            ut::KeyType key; int flavor; int state;
            bool operator==(Entry const& rhs) { return key == rhs.key;};
        };
        
        Simple const& func_;
        
        double w_, wBackup_;
        std::vector<Entry> modOps_, ops_;
        
        mutable std::vector<std::unique_ptr<double>> qkinks_;
    };
}

#include "Dynamic.impl.h"

#endif
