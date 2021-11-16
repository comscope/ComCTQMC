#ifndef INCLUDE_MPI_IOREDUCTIONS_H
#define INCLUDE_MPI_IOREDUCTIONS_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Basic.h"
#include "../io/Tensor.h"

namespace mpi {

    template <typename T>
    void pack_unpack_function_tensor(std::vector<io::Tensor<T>> & tensor, std::vector<T> & packed_vector, std::vector<std::vector<int>> const& ijkl, int const start,  bool unpack);

    template <typename T>
    void reduce_function_tensor(std::vector<io::Tensor<T>> & tensor, int const start, int root);

    template <typename T>
    void all_reduce_function_tensor(std::vector<io::Tensor<T>> & tensor, int const start);

}

#include "IOReductions.impl.h"

#endif
