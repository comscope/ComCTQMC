#ifndef INCLUDE_MEASUREMENTS_MEASUREMENTS_H
#define INCLUDE_MEASUREMENTS_MEASUREMENTS_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../io/Vector.h"
#include "../../ctqmc/include/config/Worms.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {
    
    struct All {}; struct Jackknife {};
    
    
    struct Fix {};  struct Var {};

    template<typename T, typename M>
    struct Vector {
        inline static std::string name() { return name(T(), M());};
        
        Vector() = default;
        Vector(Vector const&) = default;
        Vector(Vector&&) = default;
        Vector& operator=(Vector const&) = default;
        Vector& operator=(Vector&&) = default;
        ~Vector() = default;
        
        void add(std::vector<T> const& val, std::int64_t samples);
        
        jsx::value reduce(double fact, All, bool b64) const;
        
        jsx::value reduce(double fact, Jackknife, bool b64) const;
        
        void write(jsx::value& dest) const;
        
        std::size_t size() const;
        T at(int const i) const;
        
    private:
        
        std::int64_t samples_ = 0;
        io::Vector<T> data_;
        
        static std::string name(double const&, Fix)               { return "meas::rvecfix"; };
        static std::string name(std::complex<double> const&, Fix) { return "meas::cvecfix"; };
        
        static std::string name(double const&, Var)               { return "meas::rvecvar"; };
        static std::string name(std::complex<double> const&, Var) { return "meas::cvecvar"; };
        
        static void resize_add(io::Vector<T> const& val, io::Vector<T>& data, Fix);
        static void resize_add(io::Vector<T> const& val, io::Vector<T>& data, Var);
        
        static void resize_reduce(io::Vector<T>& data, Fix);
        static void resize_reduce(io::Vector<T>& data, Var);
    };

    
    template<typename T, typename M>
    struct Sample {
        T const& value; std::int64_t samples;
    };
    
    
    template<typename T> inline Sample<T, Fix> fix(T const& value, std::int64_t samples) { return {value, samples};};
    template<typename T> inline Sample<T, Var> var(T const& value, std::int64_t samples) { return {value, samples};};
    
    
    template<typename T, typename M>
    inline void operator<<(jsx::value& lhs, Sample<T, M> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = Vector<T, M>();
        lhs.at<Vector<T, M>>().add(std::vector<T>(1, rhs.value), rhs.samples);
    }
    
    template<typename T, typename M>
    inline void operator<<(jsx::value& lhs, Sample<std::vector<T>, M> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = Vector<T, M>();
        lhs.at<Vector<T, M>>().add(rhs.value, rhs.samples);
    }

    
    using rvecfix = Vector<double, Fix>;  using cvecfix = Vector<std::complex<double>, Fix>;
    using rvecvar = Vector<double, Var>;  using cvecvar = Vector<std::complex<double>, Var>;
    
    
    //--------------------------------------------------------------------------------------------------------------------------------
    
    
    std::int64_t reduce_steps(std::int64_t steps, All);
    std::int64_t reduce_steps(std::int64_t steps, Jackknife);
    
    
    template<typename E>
    inline void reduce(jsx::value& jOut, double fact, jsx::value const& jIn, E, bool b64) {
        if(jIn.is<rvecfix>())
            jOut = jIn.at<rvecfix>().reduce(fact, E(), b64);
        else if(jIn.is<cvecfix>())
            jOut = jIn.at<cvecfix>().reduce(fact, E(), b64);
        else if(jIn.is<rvecvar>())
            jOut = jIn.at<rvecvar>().reduce(fact, E(), b64);
        else if(jIn.is<cvecvar>())
            jOut = jIn.at<cvecvar>().reduce(fact, E(), b64);
        else if(jIn.is<jsx::object_t>()) {
            for(auto& jEntry : jIn.object()) reduce(jOut[jEntry.first], fact, jEntry.second, E(), b64);
        } else if(jIn.is<jsx::array_t>()) {
            if(!(jOut.is<jsx::array_t>() && jOut.size() == jIn.size())) jOut = jsx::array_t(jIn.size());
            int index = 0; for(auto& jEntry : jIn.array()) reduce(jOut[index++], fact, jEntry, E(), b64);
        } else
            jOut = jIn;
    }
    
    
    template<typename E>
    inline void reduce(jsx::value& jOut, jsx::value const& jIn, jsx::value const& jEtas, E, bool b64) {
        auto const pName = cfg::partition::Worm::name();
        
        jsx::value const jSign = jIn(pName)("sign").at<rvecfix>().reduce(1., E(), b64);
        auto const pSteps = reduce_steps(jIn(pName)("steps").int64(), E());
        auto const signxZp = (jSign.is<jsx::null_t>() ? 1. : jsx::at<io::rvec>(jSign).at(0))*pSteps/jEtas(pName).real64();
        
        if(!jOut.is<jsx::object_t>()) jOut = jsx::object_t();
        
        for(auto& jWorm : jIn.object()) {
            auto const wSteps = reduce_steps(jWorm.second("steps").int64(), E());
            auto const Zw = wSteps/jEtas(jWorm.first).real64();
            
            reduce(jOut[jWorm.first], Zw/signxZp, jWorm.second, E(), b64);
            
            jOut[jWorm.first]["steps"] = wSteps;
        }
        
        jOut[pName]["sign"] = jSign;
    }

    std::vector<std::string> split_by_char(std::string const& string, char const c);
    
    template  <typename measure_type, typename Value>
    void check_missing_tensor_elements(std::string const name, jsx::value& jIn);
    

    //Each worker is not guarenteed to have visited the same subspaces of the larger  worm space
    //This is not only concerning from a convergence standpoint, but it also will cause MPI to hang
    //When we try to reduce the results. So, we add a measurement of zero when there is such a miss
    template <typename Value>
    void check_missing_tensor_elements(jsx::value const& jParams, jsx::value& jMeasurements);

    void restart(jsx::value const& jIn, jsx::value & jMeas, double const eta, int const samples);

    void restart(jsx::value & jParams, jsx::value & jMeasurements);

}

#include "Measurements.impl.h"

#endif
