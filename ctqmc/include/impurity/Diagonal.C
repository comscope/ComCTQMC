
#include "Diagonal.h"
#include "../Algebra.h"

namespace imp {
    
    template<typename Mode>
    EigenValues<Mode>::EigenValues(jsx::value const& jParams, jsx::value jEigenValues, std::vector<double> const& filling, imp::Simple const* dynFunc) :
    sectorNumber_(jEigenValues.size()),
    energies_(static_cast<Energies<Mode>*>(::operator new(sizeof(Energies<Mode>)*(sectorNumber_ + 1)))) {
        mpi::cout << "Reading eigenvalues ... " << std::flush;
        
        auto const mu = jParams("mu").real64();
        
        int sector = 1;
        for(auto& jEnergies : jEigenValues.array()) {
            for(auto& energy : jsx::at<io::rvec>(jEnergies))
                energy += -mu*filling.at(sector) + (dynFunc != nullptr ? dynFunc->shift(sector) : .0);
            new(energies_ + sector++) Energies<Mode>(jParams, jsx::at<io::rvec>(jEnergies));
        }
        
        mpi::cout << "Ok" << std::endl;
    }
    
    template<typename Mode>
    EigenValues<Mode>::~EigenValues() {
        for(int s = sectorNumber_; s; --s) energies_[s].~Energies();
        ::operator delete(energies_);
    };
    
    template struct EigenValues<Host>;
#ifdef MAKE_GPU_ENABLED
    template struct EigenValues<Device>;
#endif
    
	//----------------------------------------------------------------PROPAGATOR--------------------------------------------------------------------------
    template<typename Mode>
    Propagator<Mode>::Propagator(double time, EigenValues<Mode> const& eig) : eig_(eig), time_(time), isProp_(eig_.sectorNumber() + 1), prop_(static_cast<Vector<Mode>*>(::operator new(sizeof(Vector<Mode>)*(eig_.sectorNumber() + 1)))) {}; /////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    template<typename Mode>
    Vector<Mode> const& Propagator<Mode>::at(int s) {
        if(isProp_[s]) return prop_[s];
        new(prop_ + s) Vector<Mode>(time_, eig_.at(s)); isProp_.set(s); // Soetti ok sii so, oder ???
        return prop_[s];
    };
    
    template<typename Mode>
    Vector<Mode> const& Propagator<Mode>::at(int s) const {
        if(!isProp_[s]) throw std::runtime_error("imp::Propagator::at: null pointer");
        return prop_[s];
    };
    
    
    template<typename Mode>
    void Propagator<Mode>::add(SectorNormPtrs& norms) const {
        for(SectorNormPtrs::iterator it = norms.begin; it != norms.end; ++it)
            (*it)->norm += time_*eig_.at((*it)->sector).min();
    };
    
    template<typename Mode>
    Propagator<Mode>::~Propagator() {
        if(isProp_.any()) for(int s = eig_.sectorNumber(); s; --s) if(isProp_[s]) prop_[s].~Vector();
        ::operator delete(prop_);
    };
    
    template struct Propagator<Host>;
#ifdef MAKE_GPU_ENABLED
    template struct Propagator<Device>;
#endif
    
}

