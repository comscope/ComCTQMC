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
            
            std::vector<io::cmat> green_function_on_full_axis(std::vector<io::cmat>const& green){
                std::vector<io::cmat> r(2*green.size(),io::cmat(green[0].I(),green[0].J()));
                
                auto it_forward = r.begin() + green.size();
                auto it_backward = r.begin() + green.size()-1;
                for(auto const& x : green){
                    *it_forward++ = x;
                    *it_backward-- = x.conj();
                }
                    
                return r;
            }
            
        }
        
    }
    
}

