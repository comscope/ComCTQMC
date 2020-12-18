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
