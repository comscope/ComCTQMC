#include "Susc.h"

namespace evalsim {
    
    namespace worm {
        
        namespace func {
            
            namespace susc {
                
                namespace ph {
                
                    template<typename Value>
                    jsx::value non_impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables){
                        
                        func::BosonFrequencies<Value> frequencies(jWorm);
                        
                        ////Hybridization function
                        
                        mpi::cout << "Reading hybridisation function ... " << std::flush;
                        
                        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                        std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
                        
                        std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
                        
                        mpi::cout << "Ok" << std::endl;
                                              
                                              
                        mpi::cout << "Reading function ... " << std::flush;
                        
                        std::vector<io::ctens> susc = meas::read_tensor_functions<Value,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
                        
                        mpi::cout << "Ok" << std::endl;
                        
                        std::vector<io::ctens> susc_symm(susc.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
                        func::susc::ph::enforce_symmetries<Value>(jParams,jWorm,susc,susc_symm);
                        
                        jsx::value jObservablesOut;
                        jObservablesOut["greens function"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
                        
                        mpi::cout << "Computing susceptibility ... " << std::flush;
                        
                        func::susc::ph::compute_and_subtract_disconnected<Value>(jParams,jObservables("partition")("occupation"), frequencies, susc);
                        
                        mpi::cout << "Ok" << std::endl;
                        
                        mpi::cout << "Enforcing symmetries ... " << std::flush;
                        
                        func::susc::ph::enforce_symmetries<Value>(jParams,jWorm,susc,susc_symm);
                        
                        mpi::cout << "Ok" << std::endl;
                        
                        jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
                        
                        return jObservablesOut;
                    }
                
                    template<typename Value>
                    jsx::value impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables){
                        
                        throw std::runtime_error("improved estimators are not implemented for susc worms\n");
                        
                    }
                    
                    //instantiations
                    template jsx::value non_impr_est_evalsim<double>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    template jsx::value impr_est_evalsim<double>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    template jsx::value non_impr_est_evalsim<ut::complex>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    template jsx::value impr_est_evalsim<ut::complex>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    //------------//
                    
                    template <typename Value>
                    void compute_and_subtract_disconnected(jsx::value const& jParams, jsx::value const& jOccupation, BosonFrequencies<Value> const& frequencies, std::vector<io::ctens>& full_in_connected_out){
                        
                        double const beta = jParams("beta").real64();
                        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                        auto const omega = frequencies.omega_b();
                        
                        for(std::size_t n = 0; n < full_in_connected_out.size(); ++n){
                            for(auto const& ijkl : full_in_connected_out[n].ijkl()){
                            
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                                
                                full_in_connected_out[n](i,j,k,l) *= -1;
                                                
                                if (!omega(n) and i==j and k==l){
                                    std::string const entry_ij = jHybMatrix(i)(j).string();
                                    std::string const entry_kl = jHybMatrix(k)(l).string();
                                    auto const& occ_ij = jsx::at<io::Vector<Value>>(jOccupation(entry_ij));
                                    auto const& occ_kl = jsx::at<io::Vector<Value>>(jOccupation(entry_kl));
                                    full_in_connected_out[n](i,j,k,l) -= beta*(occ_ij[0])*(occ_kl[0]);
                                }
                            }
                        }
                    }
                    
                    template <typename Value>
                    void enforce_symmetries(jsx::value const& jParams, jsx::value const& jWorm, std::vector<io::ctens> const& no_symm, std::vector<io::ctens>& symm){
                        
                        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                        
                        int const pos_and_neg_freq = (std::is_same<Value,double>::value ? 1 : 2);
                        int const nMatGB = pos_and_neg_freq*jWorm("cutoff").int64() - pos_and_neg_freq +1;
                        func::OmegaMap omega_b(nMatGB,true,!std::is_same<Value,double>::value);
                        
                        for(std::size_t n = 0; n < no_symm.size(); ++n){
                            int const m = std::is_same<Value,double>::value ? n : omega_b.pos(-omega_b(n));
                            
                            for(auto const& ijkl : no_symm[n].ijkl()){
                            
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                            
                                //Original, two fermionic swaps, hermetian symm
                                symm[n].emplace(i,j,k,l,
                                                no_symm[n].entry(i,j,k,l),
                                                0.25*(no_symm[n](i,j,k,l) + (no_symm[m](l,k,j,i)) + std::conj(no_symm[n](j,i,l,k)) + std::conj(no_symm[m](k,l,i,j)))
                                                );
                                            
                                }
                        }
                        
                    }
                    
                }
                
                namespace pp {
                
                    template<typename Value>
                    jsx::value non_impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables){
                        
                        ////Hybridization function
                        
                        mpi::cout << "Reading hybridisation function ... " << std::flush;
                        
                        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                        std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
                        
                        std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
                        
                        mpi::cout << "Ok" << std::endl;
                        
                        
                        mpi::cout << "Reading function ... " << std::flush;
                        
                        std::vector<io::ctens> susc = meas::read_tensor_functions<Value,Boson>(jMeasurements, jParams, jWorm, jHybMatrix, hyb.size());
                        
                        mpi::cout << "Ok" << std::endl;
                        
                        // Particle particle susceptibility is equal to the two-particle green function
                        mpi::cout << "Enforcing symmetries ... " << std::flush;
                        
                        std::vector<io::ctens> susc_symm(susc.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
                        func::susc::pp::enforce_symmetries<Value>(jParams,susc,susc_symm);
                    
                        mpi::cout << "Ok" << std::endl;
                        
                        jsx::value jObservablesOut;
                        
                        jObservablesOut["susceptibility"] = func::write_functions<Value>(jParams, jHybMatrix, susc_symm);
                        
                        return jObservablesOut;
                        
                    }
                
                    template<typename Value>
                    jsx::value impr_est_evalsim(jsx::value jParams, jsx::value const& jWorm, jsx::value const& jMeasurements, jsx::value const& jPartition, jsx::value const& jObservables){
                        
                        throw std::runtime_error("improved estimators are not implemented for susc worms\n");
                    }
                    
                    //instantiations
                    template jsx::value non_impr_est_evalsim<double>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    template jsx::value impr_est_evalsim<double>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    template jsx::value non_impr_est_evalsim<ut::complex>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    template jsx::value impr_est_evalsim<ut::complex>(jsx::value, jsx::value const&, jsx::value const&, jsx::value const&, jsx::value const&);
                    //------------//
                    
                    template <typename Value>
                    void compute_and_subtract_disconnected(jsx::value const& jParams, jsx::value const& jOccupation, std::vector<io::ctens>& full_in_connected_out){ std::cout << "No disconnected part of susc_pp!" << std::endl; }
                    
                    template <typename Value>
                    void enforce_symmetries(jsx::value const& jParams, std::vector<io::ctens> const& no_symm, std::vector<io::ctens>& symm){
                        
                        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                        
                        for(std::size_t n = 0; n < no_symm.size(); ++n)
                            for(auto const& ijkl : no_symm[n].ijkl()){
                        
                                auto const i = ijkl[1];
                                auto const j = ijkl[2];
                                auto const k = ijkl[3];
                                auto const l = ijkl[4];
                                            
                                symm[n].emplace(i,j,k,l,
                                                no_symm[n].entry(i,j,k,l),
                                                0.25*(no_symm[n](i,j,k,l) - no_symm[n](i,j,l,k) - no_symm[n](j,i,k,l) + no_symm[n](j,i,l,k))
                                                );
                            }
                    }
                    
                }
                
            }
            
        }
        
    }
    
}




