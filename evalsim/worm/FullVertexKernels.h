#ifndef EVALSIM_WORM_FULL_VERTEX_KERNELS
#define EVALSIM_WORM_FULL_VERTEX_KERNELS

#include "../../include/JsonX.h"
#include "../../ctqmc/include/config/Worms.h"
#include "../../ctqmc/include/Params.h"
#include "../../ctqmc/include/impurity/Tensor.h"

#include "include/Common.h"
#include "include/functions/Functions.h"
#include "include/functions/Measurements.h"
#include "include/functions/Utilities.h"

#include "../partition/Evalsim.h"

namespace evalsim {
    
    namespace worm {
        
        struct Kernels{
            static std::string const name;
        };
        
        template<typename Value>
        struct InteractionTensor {
            InteractionTensor() = delete;
            InteractionTensor(imp::Tensor<Value> const& tensor, int N) :
            N_(N),
            tensor_(N*N*N*N) {
                
                for (int i=0; i<N; i++)
                for (int j=0; j<N; j++)
                for (int k=0; k<N; k++)
                for (int l=0; l<N; l++)
                    tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l] = tensor(i,j,l,k);
                
            };
            
            Value operator()(int i, int j, int k, int l) const {
                return tensor_[N_*N_*N_*i + N_*N_*j + N_*k + l];
            };
            
        private:
            int const N_;
            io::Vector<Value> tensor_;
            
        };
        
        //We must swap operator orders to get everything into the same order as the full vertex
        //The following options do exactly that
        std::vector<io::ctens> rearrange_susc_ph(std::vector<io::ctens> const& susc_ph);
    
        std::vector<io::ctens> rearrange_susc_tph(std::vector<io::ctens> const& susc_ph);
    
        std::vector<io::ctens> rearrange_susc_pp(std::vector<io::ctens> const& susc_pp);
    
        std::vector<io::ctens> rearrange_hedin_ph(std::vector<io::ctens> const& hedin_ph);
    
        std::vector<io::ctens> rearrange_hedin_tph(std::vector<io::ctens> const& hedin_ph);
    
        std::vector<io::ctens> rearrange_hedin_pp(std::vector<io::ctens> const& hedin_pp);
        
        template<typename Value>
        jsx::value evaluateKernels(jsx::value jParams, jsx::value const& jObservables);
        
        
        template<typename Value>
        jsx::value evaluateFullVertexFromKernels(jsx::value jParams, jsx::value const& jObservables);
    
        std::vector<std::vector<int>> construct_ijkls(int const);
        
        
    }
    
}


#endif









