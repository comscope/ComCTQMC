#ifndef INCLUDE_IO_VECTOR_H
#define INCLUDE_IO_VECTOR_H

#include <vector>
#include <complex>

#include "Base64.h"
#include "../JsonX.h"

//scheiss b64 member variable !! Kack loesig im moment

namespace io {
    
template<typename T, typename std::enable_if<!std::is_same<std::complex<double>, T>::value, int>::type = 0>   // sfinae so that it is taken for ivec and rvec. However, not sure anymore why function overload below is not sufficient
jsx::value encode(std::vector<T> const& source, bool b64);

template<typename T, typename std::enable_if<!std::is_same<std::complex<double>, T>::value, int>::type = 0>
void decode(jsx::value const& source, std::vector<T>& dest);


jsx::value encode(std::vector<std::complex<double>> const& source, bool b64);
void decode(jsx::value const& source, std::vector<std::complex<double>>& dest);

    
    template<typename T> struct Vector : std::vector<T> {
        inline static std::string name() { return name(T());};
        
        Vector() = default;
        Vector(Vector const&) = default;
        Vector(Vector&&) = default;
        Vector& operator=(Vector const&) = default;
        Vector& operator=(Vector&&) = default;
        ~Vector() = default;
        
        template<typename Arg, typename... Args, typename std::enable_if< !std::is_same<typename std::decay<T>::type, Vector>::value, int>::type = 0>
        Vector(Arg&& arg, Args&&... args) : std::vector<T>(std::forward<Arg>(arg), std::forward<Args>(args)...) {}
        
        template<typename Arg, typename... Args>
        Vector(std::initializer_list<Arg>&& list, Args&&... args) : std::vector<T>(std::forward<std::initializer_list<Arg>>(list), std::forward<Args>(args)...) {}
        
        template<typename Arg, typename std::enable_if< !std::is_same<typename std::decay<Arg>::type, Vector>::value, int>::type = 0>
        Vector& operator=(Arg&& arg) { std::vector<T>::operator=(std::forward<Arg>(arg)); return *this;}
        
        bool& b64() const { return b64_;};
        void read(jsx::value const& source) { decode(source, *this);};
        void write(jsx::value& dest) const { dest = encode(*this, b64_);};
    private:
        mutable bool b64_ = false;
        
        inline static std::string name(jsx::int64_t const&) { return "io::ivec";};
        inline static std::string name(double const&) { return "io::rvec";};
        inline static std::string name(std::complex<double> const&) { return "io::cvec";};
    };
    

    typedef Vector<jsx::int64_t> ivec;
    typedef Vector<double> rvec;
    typedef Vector<std::complex<double>> cvec;

};

#include "Vector.impl.h"

#endif
