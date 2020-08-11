#ifndef CTQMC_INCLUDE_IMPURITY_UTILITIES_H
#define CTQMC_INCLUDE_IMPURITY_UTILITIES_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <new>
#include <memory>

#include "../../../include/linalg/LinAlg.h"
#include "../../../include/linalg/Operators.h"
#include "../../../include/mpi/Utilities.h"


namespace imp {
    
    template <typename Value>
    io::rvec blockNorms(jsx::value const& jOperator){
        
        io::rvec norms(jOperator.size());
        
        std::size_t index = 0;
        for(auto& jBloc : jOperator.array()) {
            if(!jBloc("target").is<jsx::null_t>()) {
                
                auto const& matrix = jsx::at<io::Matrix<Value>>(jBloc("matrix"));
                double const norm = linalg::spectral_norm(matrix);
                norms[index] = norm;
                
            }
            index++;
        }
        
        return norms;
    }

    struct NormCollection{
        
        NormCollection() = delete;
        NormCollection(std::size_t const N) : norms_(N), normsDagg_(N){}
        NormCollection(std::vector<io::rvec> norms, std::vector<io::rvec> normsDagg) : norms_(norms), normsDagg_(normsDagg){}
        
        std::vector<io::rvec> const& norms() const {return norms_;}
        std::vector<io::rvec> const& normsDagg() const {return normsDagg_;}
        
    private:
        std::vector<io::rvec> norms_;
        std::vector<io::rvec> normsDagg_;
        
    };

    template <typename Value>
    NormCollection gatherNorms(jsx::value const& jMPI, jsx::value const& jOperators){
        
        //computational bottleneck is the computation of spectral norms of operators.
        //This is often a memory bound operation -- it is better if it is split across nodes rather than ranks
        std::size_t const N = jOperators.size();
        std::size_t const size = N;
        std::size_t const rankOnNode = jMPI("rank on node").int64();
        std::size_t const nodeRank = jMPI("rank of node").int64();
        std::size_t const nNodes = jMPI("number of nodes").int64();
        std::size_t const chunk = (size + nNodes - 1)/nNodes;
        
        
        if(nNodes){//don't bother if there is only one node -- tends to just slow down computation
            auto const& jOps = jOperators.array();
            std::vector<int> sizes(N);
            std::vector<int> sizesDagg(N);
            std::vector<io::rvec> norms(N);
            std::vector<io::rvec> normsDagg(N);
            
            if (!rankOnNode) //only the head rank on each node works
                for(std::size_t index = nodeRank*chunk; index < chunk*(nodeRank + 1); ++index){
                    if (index>=size) break;
                    
                    auto const& jOp = jOps[index];
                    jsx::value const jOpDagg = linalg::conj<Value>(jOp);
                    
                    sizes[index] = jOp.size();
                    norms[index] = blockNorms<Value>(jOp);
                    
                    sizesDagg[index] = jOpDagg.size();
                    normsDagg[index] = blockNorms<Value>(jOpDagg);
                }
            
            mpi::barrier(); //extra barriers to deal with some mpi issues occasionally cropping up on systems
            mpi::all_reduce<mpi::op::sum>(sizes);
            mpi::all_reduce<mpi::op::sum>(sizesDagg);
            
            for (std::size_t i=0; i<N; i++){
                
                if(!norms[i].size()) norms[i].resize(sizes[i]);
                mpi::barrier();
                mpi::all_reduce<mpi::op::sum>(norms[i]);
                
                if(!normsDagg[i].size()) normsDagg[i].resize(sizesDagg[i]);
                mpi::barrier();
                mpi::all_reduce<mpi::op::sum>(normsDagg[i]);
                
            }
            
            return NormCollection(norms,normsDagg);
        }
        
        return NormCollection(N); // empty collection
        
    }
}

#endif  
