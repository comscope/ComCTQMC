#include "Vector.h"

namespace io {
    
    jsx::value encode(std::vector<std::complex<double>> const& source, bool b64) {
        jsx::value jDest;
        std::vector<double> real, imag;
        
        for(auto const& z : source) {
            real.push_back(z.real()); imag.push_back(z.imag());
        };
        
        jDest["real"] = encode(real, b64); jDest["imag"] = encode(imag, b64);
        
        return jDest;
    };

    
    void decode(jsx::value const& source, std::vector<std::complex<double>>& dest) {
        std::vector<double> real; decode(source("real"), real);
        std::vector<double> imag; decode(source("imag"), imag);
        
        if(real.size() != imag.size()) throw std::runtime_error("io::decode: invalid format");
        
        for(std::size_t n = 0; n < real.size(); ++n)
            dest.push_back({real[n], imag[n]});
    };

};
