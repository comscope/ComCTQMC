#ifndef INCLUDE_MPI_BASIC_IMPL_H
#define INCLUDE_MPI_BASIC_IMPL_H

#include "Basic.h"

//--------------------------------------------------schö tö tiä tü mö tiä par la barbischätöööötötötötö-------------------------------------------------------

namespace mpi {
    
    inline int rank() {
        int temp = 0;
#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &temp);
#endif
        return temp;
    };
    
    inline int number_of_workers() {
        int temp = 1;
#ifdef HAVE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &temp);
#endif
        return temp;
    };
    
    inline int processor_name_size() {
#ifdef HAVE_MPI
        return MPI_MAX_PROCESSOR_NAME;
#else
        return 1;
#endif
    }

    inline std::vector<char> processor_name() {
#ifdef HAVE_MPI 
        std::vector<char> name(processor_name_size(), '\0'); int length;
        MPI_Get_processor_name(name.data(), &length);
        return name;
#else
        return { '\0' };
#endif
    }
    
    inline void barrier() {
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    };
    
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type>
    void reduce(T& arg, int root) {
#ifdef HAVE_MPI
        T result; MPI_Reduce(&arg, &result, 1, get_data_type(T()), get_op(Op()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = result;
#endif
    };
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type>
    void reduce(std::vector<T>& arg, int root) {
#ifdef HAVE_MPI
        std::vector<T> result(rank() == root ? arg.size() : 0);
        MPI_Reduce(arg.data(), result.data(), arg.size(), get_data_type(T()), get_op(Op()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = std::move(result);
#endif
    };
    
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type>
    void all_reduce(T& arg) {
#ifdef HAVE_MPI
        T result; MPI_Allreduce(&arg, &result, 1, get_data_type(T()), get_op(Op()), MPI_COMM_WORLD);
        arg = result;
#endif
    };
    
    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type>
    void all_reduce(std::vector<T>& arg) {
#ifdef HAVE_MPI
        std::vector<T> result(arg.size());
        MPI_Allreduce(arg.data(), result.data(), arg.size(), get_data_type(T()), get_op(Op()), MPI_COMM_WORLD);
        arg = std::move(result);
#endif
    };
 
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void bcast(T& arg, int root) {
#ifdef HAVE_MPI
        MPI_Bcast(&arg, 1, get_data_type(T()), root, MPI_COMM_WORLD);
#endif
    };
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void bcast(std::vector<T>& arg, int root) {
#ifdef HAVE_MPI
        MPI_Bcast(arg.data(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
#endif
    };
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void bcast(std::basic_string<T>& arg, int root) {
#ifdef HAVE_MPI
        MPI_Bcast(&arg.front(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
#endif
    };
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void bcast(std::vector<std::basic_string<T>>& arg, int root) {
#ifdef HAVE_MPI
        for (auto & str : arg) bcast(str,root);
#endif
    };


    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void gather(std::vector<T>& arg, int root) {
#ifdef HAVE_MPI
        std::vector<T> result(rank() == root ? arg.size()*mpi::number_of_workers() : 0);
        MPI_Gather(arg.data(), arg.size(), get_data_type(T()), result.data(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = std::move(result);
#endif
    };


template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void gather(std::basic_string<T>& arg, int root) {
#ifdef HAVE_MPI
        //should we check that strings are of the same size?
        std::basic_string<T> result( rank() == root ? arg.size()*mpi::number_of_workers() : 0, ' ');
        MPI_Gather(&arg.front(), arg.size(), get_data_type(T()), &result.front(), arg.size(), get_data_type(T()), root, MPI_COMM_WORLD);
        if(rank() == root) arg = std::move(result);
#endif
    };
    
    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type>
    void scatter(std::vector<T>& arg, int size, int root) {
        if(rank() == root && size*number_of_workers() != static_cast<int>(arg.size()))
            throw std::runtime_error("mpi::scatter: send buffer has wrong size");        
#ifdef HAVE_MPI
        std::vector<T> result(size);
        MPI_Scatter(arg.data(), size, get_data_type(T()), result.data(), size, get_data_type(T()), root, MPI_COMM_WORLD);
        arg = std::move(result);
#endif
    };
    
}


#endif
