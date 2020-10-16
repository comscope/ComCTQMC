#ifndef INCLUDE_MPI_BASIC_H
#define INCLUDE_MPI_BASIC_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <string>
#include <vector>
#include <complex>

//--------------------------------------------------schö tö tiä tü mö tiä par la barbischätöööötötötötö-------------------------------------------------------

namespace mpi {
    
    enum : unsigned {
        empty       = 0x0,
        fundamental = 0x1,
        arithmetic  = 0x2,
        ordered     = 0x4,
        logical     = 0x8,
        bitwise     = 0x10
    };
    
    
    template<unsigned p>
    struct has_property {
        enum : unsigned { value = p };
    };
    
    
    template<typename T> struct property              :  has_property< empty > {};

    template<> struct property<char>                  :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<signed char>           :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned char>         :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<short>                 :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned short>        :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<int>                   :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned>              :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<long>                  :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned long>         :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<long long>             :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<unsigned long long>    :  has_property< fundamental | arithmetic | ordered | bitwise > {};
    template<> struct property<float>                 :  has_property< fundamental | arithmetic | ordered > {};
    template<> struct property<double>                :  has_property< fundamental | arithmetic | ordered > {};
    template<> struct property<std::complex<float>>   :  has_property< fundamental | arithmetic > {};
    template<> struct property<std::complex<double>>  :  has_property< fundamental | arithmetic > {};
    
    
    template<typename T, unsigned m>
    struct data_compatible_if {
        enum : bool { value = ( (property<T>::value & m) == m ) };
    };
    
    
    namespace op {
        struct min {};  struct max {}; struct sum {};  struct prod {};
        struct land {}; struct lor {}; struct lxor {};
        struct band {}; struct bor {}; struct bxor {};
    };
    
    
    template<typename T, typename Op > struct data_op_compatible { enum : bool { value = false }; };
    
    template<typename T> struct data_op_compatible<T, op::min>   :  data_compatible_if< T, fundamental | ordered > {};
    template<typename T> struct data_op_compatible<T, op::max>   :  data_compatible_if< T, fundamental | ordered > {};
    template<typename T> struct data_op_compatible<T, op::sum>   :  data_compatible_if< T, fundamental | arithmetic > {};
    template<typename T> struct data_op_compatible<T, op::prod>  :  data_compatible_if< T, fundamental | arithmetic > {};
    template<typename T> struct data_op_compatible<T, op::land>  :  data_compatible_if< T, fundamental | logical > {};
    template<typename T> struct data_op_compatible<T, op::lor>   :  data_compatible_if< T, fundamental | logical > {};
    template<typename T> struct data_op_compatible<T, op::lxor>  :  data_compatible_if< T, fundamental | logical > {};
    template<typename T> struct data_op_compatible<T, op::band>  :  data_compatible_if< T, fundamental | bitwise > {};
    template<typename T> struct data_op_compatible<T, op::bor>   :  data_compatible_if< T, fundamental | bitwise > {};
    template<typename T> struct data_op_compatible<T, op::bxor>  :  data_compatible_if< T, fundamental | bitwise > {};
    
    
#ifdef HAVE_MPI
    inline MPI_Op get_op(op::min const&)   { return MPI_MIN; };
    inline MPI_Op get_op(op::max const&)   { return MPI_MAX; };
    inline MPI_Op get_op(op::sum const&)   { return MPI_SUM; };
    inline MPI_Op get_op(op::prod const&)  { return MPI_PROD; };
    inline MPI_Op get_op(op::land const&)  { return MPI_LAND; };
    inline MPI_Op get_op(op::lor const&)   { return MPI_LOR; };
    inline MPI_Op get_op(op::lxor const&)  { return MPI_LXOR; };
    inline MPI_Op get_op(op::band const&)  { return MPI_BAND; };
    inline MPI_Op get_op(op::bor const&)   { return MPI_BOR; };
    inline MPI_Op get_op(op::bxor const&)  { return MPI_BXOR; };
    
    template<typename T> MPI_Datatype get_data_type(T const&); // redundant because of sfinae, however, sometimes redundancy is good ...
    
    inline MPI_Datatype get_data_type(char const&)                 { return MPI_CHAR; };
    inline MPI_Datatype get_data_type(signed char const&)          { return MPI_SIGNED_CHAR; };
    inline MPI_Datatype get_data_type(unsigned char const&)        { return MPI_UNSIGNED_CHAR; };
    inline MPI_Datatype get_data_type(short const&)                { return MPI_SHORT; };
    inline MPI_Datatype get_data_type(unsigned short const&)       { return MPI_UNSIGNED_SHORT; };
    inline MPI_Datatype get_data_type(int const&)                  { return MPI_INT; };
    inline MPI_Datatype get_data_type(unsigned const&)             { return MPI_UNSIGNED; };
    inline MPI_Datatype get_data_type(long const&)                 { return MPI_LONG; };
    inline MPI_Datatype get_data_type(unsigned long const&)        { return MPI_UNSIGNED_LONG; };
    inline MPI_Datatype get_data_type(long long const&)            { return MPI_LONG_LONG; };
    inline MPI_Datatype get_data_type(unsigned long long const&)   { return MPI_UNSIGNED_LONG_LONG; };
    inline MPI_Datatype get_data_type(float const&)                { return MPI_FLOAT; };
    inline MPI_Datatype get_data_type(double const&)               { return MPI_DOUBLE; };
    inline MPI_Datatype get_data_type(std::complex<float> const&)  { return MPI_C_FLOAT_COMPLEX; };
    inline MPI_Datatype get_data_type(std::complex<double> const&) { return MPI_C_DOUBLE_COMPLEX; };
#endif
    

    int const master = 0;

    inline int rank();

    inline int number_of_workers();

    inline int processor_name_size();

    inline std::vector<char> processor_name();

    inline void barrier();

    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void reduce(T& arg, int root);

    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void reduce(std::vector<T>& arg, int root);

    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void all_reduce(T& arg);

    template<typename Op, typename T, typename std::enable_if<data_op_compatible<T, Op>::value, int>::type = 0>
    void all_reduce(std::vector<T>& arg);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(T& arg, int root);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(std::vector<T>& arg, int root);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(std::basic_string<T>& arg, int root);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void bcast(std::vector<std::basic_string<T>>& arg, int root);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void gather(std::vector<T>& arg, int root);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void gather(std::basic_string<T>& arg, int root);

    template<typename T, typename std::enable_if<data_compatible_if< T, fundamental >::value, int>::type = 0>
    void scatter(std::vector<T>& arg, int size, int root);
    
}

#include "Basic.impl.h"

#endif
