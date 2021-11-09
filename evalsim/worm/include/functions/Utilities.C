#include "Utilities.h"

namespace evalsim {
    
    namespace worm {
        
        namespace func {
            
            iOmega::iOmega(double const beta, int const is_fermionic) : beta_(beta), is_fermionic_(is_fermionic) {};
            
            std::complex<double> iOmega::operator()(int n) const {
                return {.0, (2*n + is_fermionic_)*M_PI/beta_};
            };

            OmegaMap::OmegaMap(std::size_t const n, bool bosonic, bool symmetric) : n_(n), val_(n,0){
                    
                if (symmetric){
                    if (!bosonic and !n%2) throw std::runtime_error("symmetric OmegaMap: number of fermionic frequencies must be even\n");
                    if (bosonic and n%2 == 0) throw std::runtime_error("symmetric OmegaMap: number of bosonic frequencies must be odd\n");
                }
                
                int const start = !symmetric ? 0 : -n_/2;
                int const shift = bosonic ? 0 : 1;
                
                for (std::size_t i=0;i<n_; i++){
                    val_[i]=2*(start+i)+shift;
                    pos_[val_[i]]=i;
                }
                    
            }
            
            std::vector<io::cmat> half_to_full(std::vector<io::cmat> const& half){
                std::vector<io::cmat> full(2*half.size(), io::cmat(half[0].I(), half[0].J()));
                
                auto it_forward = full.begin() + half.size();
                auto it_backward = full.begin() + half.size()-1;
                
                for(int n=0; n<half.size(); n++){
                    auto const& x = half[n];
                    *it_forward++ = x;
                    *it_backward-- = x.conj(); //TODO: should be transposed too
                }
                
                return full;
            }
            
            template<typename Value>
            std::vector<io::cmat> get_green_from_obs(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name){
                
                //use the green imprsum measurement if available (assume user has a good reason).
                //Otherwise, use the partition function measurement (which is typically better).
                return jObservables.is(cfg::green_imprsum::Worm::name()) ?
                meas::read_matrix_functions_from_obs<ut::complex>(jObservables(cfg::green_imprsum::Worm::name())(func_name), jParams, jParams(cfg::green_imprsum::Worm::name()), jHybMatrix, hybSize) :
                meas::read_matrix_functions_from_obs<ut::complex>(jObservables(cfg::partition::Worm::name())(func_name), jParams, jParams(cfg::partition::Worm::name()), jHybMatrix, hybSize);
                
            }
            
            template<typename Value>
            std::vector<io::Matrix<Value>> get_moments_from_obs(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name){
                
                //use the green imprsum measurement if available (assume user has a good reason).
                //Otherwise, use the partition function measurement (which is typically better).
                return jObservables.is(cfg::green_imprsum::Worm::name()) ?
                meas::read_matrix_moments_from_obs<Value>(jObservables(cfg::green_imprsum::Worm::name())(func_name), jParams, jParams(cfg::green_imprsum::Worm::name()), jHybMatrix, hybSize) :
                meas::read_matrix_moments_from_obs<Value>(jObservables(cfg::partition::Worm::name())(func_name), jParams, jParams(cfg::partition::Worm::name()), jHybMatrix, hybSize);
                
            }
            
            template <typename Value>
            void extend_hyb(std::vector<io::cmat>& hyb, std::vector<io::Matrix<Value>> const& hybMoments, int const nf, double const beta){
                
                auto const I = hyb[0].I();
                auto const N = hyb.size();
                
                iOmega const iomega(beta);
                
                for (int n=N; n<nf; n++){
                    io::cmat tmp(I,I);
                    auto const w = iomega(n);
                    
                    for (int i=0; i<I; i++)
                        for (int j=0; j<I; j++){
                            if (std::abs(hybMoments[0](i,j)) > 0)
                                tmp(i,j) = hybMoments[0](i,j)/(w - hybMoments[1](i,j));
                            
                        }
                    
                    hyb.push_back(tmp);
                }
            }
            
            template <typename Value>
            void green_and_self_function_on_full_axis(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix,
                                                      std::vector<io::cmat>& hyb, std::vector<io::Matrix<Value>> const& hybMoments, int const nf,
                                                      std::vector<io::cmat>& green, std::vector<io::cmat>& self){
                
                if (green.size() != self.size() and nf != green.size()) throw std::runtime_error("Green's and self-energy functions are not the correct sizes");
                int const mid = int(nf/2);
                
                auto const beta = jParams("beta").real64();
                auto green_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
                auto self_half_axis = func::get_green_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"self-energy");
                
                if (mid > hyb.size()){
                    auto const green_moments = func::get_moments_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"green");
                    auto const self_moments = func::get_moments_from_obs<Value>(jParams,jObservables,jHybMatrix,hyb.size(),"self-energy");
                    
                    extend_hyb(hyb, hybMoments, mid, jParams("beta").real64());
                    partition::func::add_self_tail<Value>(jParams, self_half_axis, self_moments, mid);
                    partition::func::add_green_tail<Value>(jParams, hyb, self_half_axis, green_half_axis);
                }
                
                green = half_to_full(green_half_axis);
                self = half_to_full(self_half_axis);
                
            }
            
            template void green_and_self_function_on_full_axis<double>(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix,
                                                                       std::vector<io::cmat> & hyb, std::vector<io::Matrix<double>> const& hybMoments, int const nf,
                                                                       std::vector<io::cmat>& green, std::vector<io::cmat>& self);
            template void green_and_self_function_on_full_axis<ut::complex>(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix,
                                                                            std::vector<io::cmat> & hyb, std::vector<io::Matrix<ut::complex>> const& hybMoments, int const nf,
                                                                            std::vector<io::cmat>& green, std::vector<io::cmat>& self);
            
            template std::vector<io::cmat> get_green_from_obs<double>(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name);
            template std::vector<io::cmat> get_green_from_obs<ut::complex>(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name);
            
            
            template std::vector<io::Matrix<double>> get_moments_from_obs<double>(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name);
            template std::vector<io::Matrix<ut::complex>> get_moments_from_obs<ut::complex>(jsx::value const& jParams, jsx::value const& jObservables, jsx::value const& jHybMatrix, std::size_t hybSize, std::string const func_name);
            
            
            template void extend_hyb<double>(std::vector<io::cmat>& hyb, std::vector<io::rmat> const& hybMoments, int const nf, double const beta);
            template void extend_hyb<ut::complex>(std::vector<io::cmat>& hyb, std::vector<io::cmat> const& hybMoments, int const nf, double const beta);
            
            
        }
        
    }
    
}

