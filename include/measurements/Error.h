#ifndef INCLUDE_MEASUREMENTS_ERROR_H
#define INCLUDE_MEASUREMENTS_ERROR_H

#include <vector>
#include <complex>
#include <cstring>

#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../io/Vector.h"

#include "Measurements.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {


    void error(std::vector<double>& arg, Jackknife);
    void error(std::vector<double>& arg, Variance);

    template<typename E>
    inline void error(jsx::value& jArg, E) {
        if(jArg.is<io::rvec>())
            error(jArg.at<io::rvec>(), E());
        else if(jArg.is<io::cvec>()) {
            auto& arg = jArg.at<io::cvec>();
            
            std::vector<double> real, imag;
            for(auto& entry : arg) {
                real.push_back(entry.real()); imag.push_back(entry.imag());
            }
            error(real, E()); error(imag, E());
            for(std::size_t i = 0; i < arg.size(); ++i)
                arg[i] = {real[i], imag[i]};
        } else if(jArg.is<jsx::object_t>()) {
            for(auto& jEntry : jArg.object()) error(jEntry.second, E());
        } else if(jArg.is<jsx::array_t>()) {
            for(auto& jEntry : jArg.array()) error(jEntry, E());
        }
    }

    inline void subtract(std::vector<double>& arg1, std::vector<double> const& arg2){
        for (int i=0; i < arg1.size(); i++)
            arg1[i] -= arg2[i];
    }

    inline void subtract(jsx::value& jArg1, jsx::value const& jArg2) {
        if(jArg1.is<io::rvec>())
            subtract(jArg1.at<io::rvec>(), jArg2.at<io::rvec>());
        else if(jArg1.is<io::cvec>()) {
            auto& arg1 = jArg1.at<io::cvec>();
            auto& arg2 = jArg2.at<io::cvec>();
            
            std::vector<double> real1, imag1, real2, imag2;
            for(auto& entry : arg1) {
                real1.push_back(entry.real()); imag1.push_back(entry.imag());
            }
            
            for(auto& entry : arg2) {
                real2.push_back(entry.real()); imag2.push_back(entry.imag());
            }
            
            subtract(real1, real2); subtract(imag1, imag2);
            
            for(std::size_t i = 0; i < arg1.size(); ++i){
                arg1[i] = {real1[i], imag1[i]};
            }
            
        } else if(jArg1.is<jsx::object_t>()) {
            for(auto& jEntry : jArg1.object()){
                auto entry = jEntry.first;
                subtract( jArg1(jEntry.first), jArg2(jEntry.first) );
            }
            
        } else if(jArg1.is<jsx::array_t>()) {
            for(int i = 0; i < jArg1.array().size(); i++ )
                subtract(jArg1.array()[i], jArg2.array()[i]);
        }
    }

    struct Error {
        Error() = default;
        ~Error() = default;
        
        void add(jsx::value const& jObservable, jsx::value const& jObservable0);
        jsx::value finalize(double norm, jsx::value const& jObservable0);
        
    private:
        jsx::value jMean_, jSquare_;
        
        static void add(jsx::value& jMean, jsx::value& jSquare, jsx::value const& jObservable, jsx::value const& jObservable0);
        
        static void finalize(jsx::value& jMean, jsx::value const& jSquare, jsx::value const& jObservable0, double const norm);
    };

}

#include "Error.impl.h"

#endif
