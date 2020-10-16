
#include "Hyb.h"

namespace bath {
    
    template<typename Value>
    Fit<Value>::Fit(double const beta, std::vector<std::complex<double>> const& hyb, std::vector<std::complex<double>> const& hybTransp) {
        if(hyb.size() != hybTransp.size()) throw std::runtime_error("bath::Fit: size if conjugate hybridisation functions not the same !");
        if(!hyb.size()) throw std::runtime_error("bath::Fit: no entries in hybridisation function !");
        
        N_ = std::max(static_cast<int>(hyb.size()*.1), 1);
        ut::complex D = .0, iwD = 0, D2 = .0, iwD2 = .0;
        for(std::size_t m = hyb.size() - N_; m < hyb.size(); ++m) {
            check(hyb[m], hybTransp[m], Value());
            
            ut::complex const iw{.0, M_PI*(2*m + 1)/beta};
            D    += hyb[m] + std::conj(hybTransp[m]);  // real if real
            iwD  += iw*(hyb[m] - std::conj(hybTransp[m])); // real if real
            D2   += std::abs(hyb[m])*std::abs(hyb[m]) + std::abs(hybTransp[m])*std::abs(hybTransp[m]); // real
            iwD2 += iw*(std::abs(hyb[m])*std::abs(hyb[m]) - std::abs(hybTransp[m])*std::abs(hybTransp[m])); // zero if real
        }
        
        moment_ = ut::to_val<Value>((iwD*D2 - D*iwD2)/(2.*N_*D2 - std::abs(D)*std::abs(D)));
        eps_ = ut::to_val<Value>((2.*N_*iwD2 - iwD*std::conj(D))/(2.*N_*D2 - std::abs(D)*std::abs(D)));
    }
    
    template<typename Value>
    void Fit<Value>::check(ut::complex pos, ut::complex neg, double) {
        if( std::abs(pos - neg) > 1.e-12*(std::abs(pos) + std::abs(neg)) )
            throw std::runtime_error("bath::Fit<double>: hybridisation function is not real !");
    }
    
    template struct Fit<double>;
    template struct Fit<ut::complex>;
    
    
    template<typename Value>
    Simple<Value>::Simple(jsx::value const& jParams, std::vector<std::complex<double>> const& hyb, std::vector<std::complex<double>> hybTransp) :
    I_(std::max(static_cast<int>((jParams.is("hybridisation factor") ? jParams("hybridisation factor").real64() : 4.)*hyb.size()), 1)),
    fact_(I_/static_cast<double>(ut::KeyMax)),
    data_(I_ + 2) {
        Fit<Value> fit(ut::beta(), hyb, hybTransp);
        
        Value fact = -fit.moment()/(1. + std::exp((ut::real(fit.eps()) < .0 ? 1. : -1.)*fit.eps()*ut::beta()));
        for(std::size_t i = 0; i < I_ + 1; ++i) {
            double const time = ut::beta()*i/static_cast<double>(I_);
            data_[i] = fact*std::exp(fit.eps()*(ut::real(fit.eps()) < .0 ? ut::beta() - time : -time));
        }
        
        for(std::size_t m = 0; m < hyb.size(); ++m) {
            ut::complex const iw{.0, M_PI*(2*m + 1)/ut::beta()};
            double const smooth = (m + fit.N() < hyb.size() ? 1. : (hyb.size() - m)/static_cast<double>(fit.N()));
            ut::complex const value = (hyb[m] - fit.moment() /(iw - fit.eps()))*smooth;
            ut::complex const valueTransp = (hybTransp[m] - ut::conj(fit.moment())/(iw - ut::conj(fit.eps())))*smooth;
            
            for(std::size_t i = 0; i < I_ + 1; ++i) {
                double const time = ut::beta()*i/static_cast<double>(I_);
                ut::complex exp{std::cos(-time*iw.imag()), std::sin(-time*iw.imag())}; // speed-up possible ....
                data_[i] += ut::to_val<Value>(exp*value + std::conj(exp*valueTransp))/ut::beta();
            }
        }
        
        data_.back() = .0;
    };
    
    template<typename Value>
    Value Simple<Value>::get(ut::KeyType key) const {
        double it = fact_*key; int i0 = static_cast<int>(it);
        return (1. - (it - i0))*data_[i0] + (it - i0)*data_[i0 + 1];
    };
    
    template struct Simple<double>;
    template struct Simple<ut::complex>;
    
    
    template<typename Value>
    Hyb<Value>::Hyb(jsx::value const& jParams, jsx::value const& jMatrix, jsx::value jFunctions) :
    flavors_(2*jParams("hybridisation")("matrix").size()),  //!!!!!!!!!! test against jAtomic !!!
    matrix_(flavors_*flavors_, nullptr),
    bath_(flavors_),
    isR_(flavors_) {
        mpi::cout << "Reading in hybridisation ... " << std::flush;

        for(auto const& row : jMatrix.array())
            if(row.size() != jMatrix.size())
                throw std::runtime_error("Hyb: hybridisation is not a matrix.");

        ga::Join labels(jMatrix.size());
        for(std::size_t i = 0; i < jMatrix.size(); ++i)
            for(std::size_t j = 0; j < jMatrix.size(); ++j)
                if(jMatrix(i)(j).string() != "") {
                    if(jMatrix(j)(i).string() == "")
                        throw std::runtime_error("Hyb: invalid hybridisation matrix.");
                    
                    std::string const entry = jMatrix(i)(j).string();
                    std::string const entryTransp = jMatrix(j)(i).string();

                    if(!data_.count(entry))
                        data_.emplace(entry, Simple<Value>(jParams, jsx::at<io::cvec>(jFunctions(entry)), jsx::at<io::cvec>(jFunctions(entryTransp))));

                    matrix_[(2*i + 1)  + flavors_*2*j] = &data_.at(entry);
                    
                    labels.join(i, j);
                }
        blocks_.resize(labels.clean());

        for(std::size_t i = 0; i < jMatrix.size(); ++i)
            for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                if(labels.label(i) == labels.label(j) && jMatrix(i)(j).string() == "")
                    throw std::runtime_error("Blocks: hybridisation matrix is not a block-diagonal matrix.");
                if(labels.label(i) != labels.label(j) && jMatrix(i)(j).string() != "")
                    throw std::runtime_error("Blocks: hybridisation matrix is not a block-diagonal matrix.");
            }

        for(std::size_t i = 0; i < jMatrix.size(); ++i) {
            blocks_[labels.label(i)].flavorsL().push_back(2*i + 1);  isR_[2*i + 1] = false;
            blocks_[labels.label(i)].flavorsR().push_back(2*i    );  isR_[2*i    ] = true;
            
            bath_[2*i + 1] = bath_[2*i] = labels.label(i);
        }

        mpi::cout << "Ok" << std::endl;
    };
    
    template<typename Value>
    Value Hyb<Value>::operator()(int flavorL, int flavorR, ut::KeyType key) const {
        return key < 0 ? -matrix_[flavorL + flavors_*flavorR]->get(key + ut::KeyMax) : matrix_[flavorL + flavors_*flavorR]->get(key);
    };
    
    template struct Hyb<double>;
    template struct Hyb<ut::complex>;
    
    
    
}

