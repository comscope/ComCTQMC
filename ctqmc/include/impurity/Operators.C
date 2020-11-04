
#include "Operators.h"
#include "../Algebra.h"


namespace imp {
    
    template<typename Mode, typename Value>
    Operator<Mode,Value>::Operator(itf::EigenValues const& eigItf) :
        eig_(get<Mode>(eigItf)),
        isMap_(eig_.sectorNumber() + 1),
        isMat_(eig_.sectorNumber() + 1),
        map_(new SectorNorm[eig_.sectorNumber() + 1]),
        mat_(static_cast<Matrix<Mode, Value>*>(::operator new(sizeof(Matrix<Mode, Value>)*(eig_.sectorNumber() + 1)))) {
        };
    
    template<typename Mode, typename Value>
    Operator<Mode,Value>::Operator(char const option, itf::EigenValues const& eigItf) : Operator(eigItf) {
        if(option == '1') {
            for(int s = eig_.sectorNumber(); s; --s) {
                set_map(s) = { s, .0 };
                mat(s, typename Matrix<Mode, Value>::Identity(eig_.at(s).dim()));
            }
        } else
            throw std::runtime_error("Tr: option in operator constructor not defined");
    };
    
    template<typename Mode, typename Value>
    Operator<Mode,Value>::Operator(jsx::value const& jOperator, itf::EigenValues const& eigItf, io::rvec const& norms) : Operator(eigItf) {
        if(static_cast<int>(jOperator.size()) != eig_.sectorNumber())
            throw(std::runtime_error("Tr: wrong number of sectors."));
        
        std::vector<int> temp(eig_.sectorNumber() + 1, 0); int start_sector = 1;
        for(auto& jBloc : jOperator.array()) {
            if(!jBloc("target").is<jsx::null_t>()) {
                
                int target_sector = jBloc("target").int64() + 1;
                
                if(target_sector < 1 || eig_.sectorNumber() < target_sector)
                    throw std::runtime_error("Tr: target sector out of range.");
                if(temp[target_sector]++)
                    throw std::runtime_error("Tr: target sector not unique.");
                
                auto const& matrix = jsx::at<io::Matrix<Value>>(jBloc("matrix"));
                
                if(matrix.I() != eig_.at(target_sector).dim0() || matrix.J() != eig_.at(start_sector).dim0())
                    throw std::runtime_error("Tr: invalid matrix dimensions");
                
                double norm;
                if (norms.size()){ norm = norms[start_sector-1]; }
                else norm = linalg::spectral_norm(matrix);
                
                if(norm != .0) {
                    set_map(start_sector) = { target_sector, std::log(norm) };
                    mat(start_sector, eig_.at(target_sector).dim(), eig_.at(start_sector).dim(), matrix);
                } else
                    set_map(start_sector) = { 0, .0 };
                
                //jBloc("matrix") = jsx::empty_t(); //TODO: Freeing jOperator matrices helps with memory noticeably? should maybe be a clean(jParams) func
            } else
                set_map(start_sector) = { 0, .0 };
            ++start_sector;
        }
    };
    
    template<typename Mode, typename Value>
    Operator<Mode,Value>::~Operator() {
        if(isMat_.any()) for(int s = eig_.sectorNumber(); s; --s) if(isMat_[s]) mat_[s].~Matrix();
        ::operator delete(mat_); delete [] map_;
    };
    
    template<typename Mode, typename Value>
    int Operator<Mode,Value>::isMap(int s) const { return isMap_[s];};
    
    template<typename Mode, typename Value>
    SectorNorm& Operator<Mode,Value>::set_map(int s) { isMap_.set(s); return map_[s];};
    
    template<typename Mode, typename Value>
    SectorNorm& Operator<Mode,Value>::map(int s) {
        if(!isMap(s)) throw std::runtime_error("imp::Operator::map: invalid sector");
        return map_[s];
    };
    
    template<typename Mode, typename Value>
    SectorNorm const& Operator<Mode,Value>:: map(int s) const {
        if(!isMap(s)) throw std::runtime_error("imp::Operator::map const: invalid sector");
        return map_[s];
    };
    
    template<typename Mode, typename Value>
    int Operator<Mode,Value>::isMat(int s) const { return isMat_[s];};
    
    template<typename Mode, typename Value>
    Matrix<Mode, Value>& Operator<Mode,Value>::mat(int s) {
        if(!isMat_[s]) throw std::runtime_error("imp::Operator::mat: invalid sector");
        return mat_[s];
    };
    
    template<typename Mode, typename Value>
    Matrix<Mode, Value> const& Operator<Mode,Value>::mat(int s) const {
        if(!isMat_[s]) throw std::runtime_error("imp::Operator::mat const: invalid sector");
        return mat_[s];
    };
    
    template<typename Mode, typename Value>
    void Operator<Mode,Value>::assign(std::vector<int> const& sectors, SectorNormPtrs& arg) {
        arg.end = arg.begin;
        for(auto sec : sectors)
            if(!isMap(sec)) {
                set_map(sec) = { sec, .0 };
                *arg.end++ = map_ + sec;
            }
    };
    
    template<typename Mode, typename Value>
    int Operator<Mode,Value>::missing(SectorNormPtrs& missing, SectorNormPtrs const& requested) {
        missing.begin = missing.end = requested.end;
        for(SectorNormPtrs::iterator it = requested.begin; it != requested.end; ++it)
            if(!isMap((*it)->sector)) {
                set_map((*it)->sector) = { (*it)->sector, .0 }; // values are set by following map
                *missing.end++ = map_ + (*it)->sector;
            }
        return missing.begin != missing.end;
    };
    
    template<typename Mode, typename Value>
    void Operator<Mode,Value>::map(SectorNormPtrs& arg) const {
        SectorNormPtrs::iterator const end = arg.end; arg.end = arg.begin;
        for(SectorNormPtrs::iterator it = arg.begin; it != end; ++it) {
            (*it)->norm += map((*it)->sector).norm;
            if(((*it)->sector = map((*it)->sector).sector))
                *arg.end++ = *it;
        }
    };
    
    template struct Operator<imp::Host,double>;
    template struct Operator<imp::Host,ut::complex>;
#ifdef MAKE_GPU_ENABLED
    template struct Operator<imp::Device,double>;
    template struct Operator<imp::Device,ut::complex>;
#endif
    
    //------------------------------------------------------------------------------------------------------------------------------
    
    template<typename Mode, typename Value>
    Operators<Mode,Value>::Operators(jsx::value const& jParams, jsx::value const& jOperators, itf::EigenValues const& eigItf) :
    flavors_(2*jOperators.size()),
    ops_(static_cast<Operator<Mode, Value>*>(::operator new(flavors_*sizeof(Operator<Mode, Value>)))) {
        mpi::cout << "Reading operators ... " << std::flush;
        
        auto norms = gatherNorms<Value>(jParams("mpi structure"), jOperators);
        
        int i = 0;
        for(auto& jOp : jOperators.array()) {
            jsx::value jOpDagg = linalg::conj<Value>(jOp);
            
            new(ops_ + 2*i    ) Operator<Mode, Value>(jOp, eigItf, norms.norms()[i]);
            new(ops_ + 2*i + 1) Operator<Mode, Value>(jOpDagg, eigItf, norms.normsDagg()[i]);
            
            ++i;
        }
        
        mpi::cout << "Ok" << std::endl;
    }
    
    template<typename Mode, typename Value>
    Operators<Mode,Value>::~Operators() {
        for(int f = 0; f < flavors_; ++f) ops_[f].~Operator();
        ::operator delete(ops_);
    };
    
    template struct Operators<imp::Host,double>;
    template struct Operators<imp::Host,ut::complex>;
#ifdef MAKE_GPU_ENABLED
    template struct Operators<imp::Device,double>;
    template struct Operators<imp::Device,ut::complex>;
#endif
    
    //------------------------------------------------------------------------------------------------------------------------------
    
    
    
}

