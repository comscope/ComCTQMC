#ifndef CTQMC_INCLUDE_IMPURITY_DIAGONAL_H
#define CTQMC_INCLUDE_IMPURITY_DIAGONAL_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>


#include "Algebra.h"
#include "BitSet.h"
#include "Dynamic.h"
#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"
#include "../../../include/JsonX.h"


namespace imp {

    namespace itf {
        
        struct EigenValues {
            virtual int sectorNumber() const = 0;
            virtual ~EigenValues() = default;
        };
        
    };

    struct SectorNorm { int sector; double norm;};
    struct SectorNormPtrs { typedef SectorNorm** iterator; SectorNorm** begin; SectorNorm** end;};
    
    // Koennte von std::vector abgeleited werden ... benötigt aber definition von move assignement (nothrow !) an verschiedenen stellen ... fuck that
    template<typename Mode>
    struct EigenValues : itf::EigenValues {
        EigenValues() = delete;
        EigenValues(jsx::value const& jParams, jsx::value jEigenValues, std::vector<double> const& filling, imp::Simple const* dynFunc);
        EigenValues(EigenValues const&) = delete;
        EigenValues(EigenValues&&) = delete;
        EigenValues& operator=(EigenValues const&) = delete;
        EigenValues& operator=(EigenValues&&) = delete;        
        ~EigenValues();
        
        int sectorNumber() const { return sectorNumber_;};
        Energies<Mode> const& at(int s) const { return energies_[s]; };

    private:
        int const sectorNumber_;
        Energies<Mode>* energies_;
    };

	//----------------------------------------------------------------PROPAGATOR--------------------------------------------------------------------------
    template<typename Mode>
	struct Propagator {
        Propagator() = delete;
        Propagator(double time, EigenValues<Mode> const& eig);
        Propagator(Propagator const&) = delete;
        Propagator(Propagator&&) = delete;
        Propagator& operator=(Propagator const&) = delete;
        Propagator& operator=(Propagator&&) = delete;
        
        double time() const { return time_;};
        Vector<Mode> const& at(int s);
        Vector<Mode> const& at(int s) const;
        void add(SectorNormPtrs& norms) const;
        ~Propagator();
        
	private:
		EigenValues<Mode> const& eig_;

		double const time_;
        BitSet isProp_;     //eleganz vo arsch vo chue aber z'schnellschte won ich bis jetzt gfunde han
		Vector<Mode>* const prop_;
	};
    
    template<typename Mode> EigenValues<Mode>& get(itf::EigenValues& eigenValuesItf) {
        return static_cast<EigenValues<Mode>&>(eigenValuesItf);
    };
    
    template<typename Mode> EigenValues<Mode> const& get(itf::EigenValues const& eigenValuesItf) {
        return static_cast<EigenValues<Mode> const&>(eigenValuesItf);
    };
	
}

#include "Diagonal.impl.h"

#endif  
