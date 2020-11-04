
#include "Observables.h"
#include "../Algebra.h"

namespace imp {
    
    
    template<typename Mode, typename Value>
    BullaOperators<Mode,Value>::BullaOperators(jsx::value const& jMPI, jsx::value const& jInteraction, jsx::value const& jOperators, itf::EigenValues const& eig) :
    flavors_(2*jOperators.size()),
    ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
        mpi::cout << "Reading bulla operators ... " << std::flush;
        
        if(static_cast<int>(jInteraction.size()) != eig.sectorNumber())
            throw(std::runtime_error("imp: wrong number of sectors in interaction."));
        

        jsx::array_t jBullas(jOperators.size());
        
        int i = 0;
        for(auto& jOp : jOperators.array()) {
            
            jsx::value jBulla;
            linalg::mult<Value>('n', 'n',  1., jOp, jInteraction, .0, jBulla);
            linalg::mult<Value>('n', 'n', -1., jInteraction, jOp, 1., jBulla);
            
            jBullas[i] = std::move(jBulla);
            ++i;
        }

        auto norms = gatherNorms<Value>(jMPI, jBullas);
            
        for (i = 0; i < jOperators.size(); ++i){
            jsx::value jBullaDagg = linalg::conj<Value>(jBullas[i]);
            
            new(ops_ + 2*i    ) Operator<Mode, Value>(jBullas[i], eig, norms.norms()[i]);
            new(ops_ + 2*i + 1) Operator<Mode, Value>(jBullaDagg, eig, norms.normsDagg()[i]);
        }
        
        mpi::cout << "Ok" << std::endl;
    };
    
    template<typename Mode, typename Value>
    BullaOperators<Mode,Value>::~BullaOperators() {
        for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
        ::operator delete(ops_);
    };
    
    
    template<typename Mode, typename Value>
    Occupation<Mode,Value>::Occupation(jsx::value const& jOperators, itf::EigenValues const& eig) :
    flavors_(jOperators.size()),
    ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
        mpi::cout << "Reading occupation ... " << std::flush;
        
        int i = 0;
        for(auto& jOp : jOperators.array()) {
            jsx::value jOcc;
            linalg::mult<Value>('c', 'n', 1., jOp, jOp, .0, jOcc);
            
            new(ops_ + i) Operator<Mode, Value>(jOcc, eig);
            
            ++i;
        }
        
        mpi::cout << "Ok" << std::endl;
    };
    
    template<typename Mode, typename Value>
    Occupation<Mode,Value>::~Occupation() {
        for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
        ::operator delete(ops_);
    };
    
    template<typename Mode, typename Value>
    BullaOccupation<Mode,Value>::BullaOccupation(jsx::value const& jParams, std::vector<double> const& filling, jsx::value jEigenValues, jsx::value const& jOperators, itf::EigenValues const& eig) :
    flavors_(jOperators.size()),
    ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
        mpi::cout << "Reading bulla occupation ... " << std::flush;
        
        auto const mu = jParams("mu").real64();
        int sector = 1;
        
        for(auto& jEnergies : jEigenValues.array())
            for(auto& energy : jsx::at<io::rvec>(jEnergies))
                energy += -mu*filling.at(sector);
        
        jsx::value jMatrixEigenValues = linalg::diag_to_operator<Value>(jEigenValues);
        
        int i = 0;
        for(auto const& jOp : jOperators.array()) {
            jsx::value jOcc;
            linalg::mult<Value>('c', 'n', 1., jOp, jOp, .0, jOcc);
            
            jsx::value jBullaOcc;
            linalg::mult<Value>('n', 'n',  1., jMatrixEigenValues, jOcc, .0, jBullaOcc);
            linalg::mult<Value>('n', 'n', -1., jOcc, jMatrixEigenValues, 1., jBullaOcc);
            
            new(ops_ + i) Operator<Mode, Value>(jBullaOcc, eig);
            
            ++i;
        }
        
        mpi::cout << "Ok" << std::endl;
    };
    
    template<typename Mode, typename Value>
    BullaOccupation<Mode,Value>::~BullaOccupation() {
        for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
        ::operator delete(ops_);
    };
    
    template struct BullaOperators<imp::Host,double>;
    template struct BullaOperators<imp::Host,ut::complex>;
    
    template struct Occupation<imp::Host,double>;
    template struct Occupation<imp::Host,ut::complex>;
    
    template struct BullaOccupation<imp::Host,double>;
    template struct BullaOccupation<imp::Host,ut::complex>;
    
#ifdef MAKE_GPU_ENABLED
    template struct BullaOperators<imp::Device,double>;
    template struct BullaOperators<imp::Device,ut::complex>;
    
    template struct Occupation<imp::Device,double>;
    template struct Occupation<imp::Device,ut::complex>;
    
    template struct BullaOccupation<imp::Device,double>;
    template struct BullaOccupation<imp::Device,ut::complex>;
#endif
    
}

