#ifndef INCLUDE_MPI_IOREDUCTIONS_IMPL_H
#define INCLUDE_MPI_IOREDUCTIONS_IMPL_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif


#include "IOReductions.h"

namespace mpi {		
		

    template <typename T>
    void pack_unpack_function_tensor(std::vector<io::Tensor<T>> & tensor, std::vector<T> & packed_vector, std::vector<std::vector<int>> const& ijkls, int const start, bool unpack){
        
        auto it = packed_vector.begin();
        
        for (int n=0; n<tensor.size(); n++)
            for (auto const& ijkl : ijkls){
                auto const i = ijkl[1];
                auto const j = ijkl[2];
                auto const k = ijkl[3];
                auto const l = ijkl[4];
                
                if (unpack)
                    *it++ = tensor[n].at(i,j,k,l);
                else
                    if (tensor[n].is(i,j,k,l))
                        tensor[n](i,j,k,l) = *it++;
                    else{
                        auto const& entry = tensor[start].entry(i,j,k,l);
                        tensor[n].emplace(i,j,k,l, entry, *it++);
                    }
                }
        
    }

    template <typename op, typename T>
    void reduce_function_tensor(std::vector<io::Tensor<T>> & tensor, int const start, int root){
        
        if (number_of_workers()>1){
            auto const& ijkls = tensor[start].ijkl();
            int N = tensor.size() * ijkls.size();
            std::vector<T> packed_vector(N, 0.);
            
            pack_unpack_function_tensor(tensor,packed_vector,ijkls,start, true);
            reduce<op,T>(packed_vector, root);
            if (rank() == root)
                pack_unpack_function_tensor(tensor,packed_vector,ijkls,start, false);
        }
        
    }


    template <typename op, typename T>
    void all_reduce_function_tensor(std::vector<io::Tensor<T>> & tensor, int const start){
        
        if (number_of_workers()>1){
            auto const& ijkls = tensor[start].ijkl();
            int N = tensor.size() * ijkls.size();
            std::vector<T> packed_vector(N, 0.);
            
            pack_unpack_function_tensor(tensor,packed_vector,ijkls,start, true);
            all_reduce<op,T>(packed_vector);
            pack_unpack_function_tensor(tensor,packed_vector,ijkls,start, false);
        }
        
    }

}

#endif
