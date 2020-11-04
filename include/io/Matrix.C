#include "Matrix.h"

namespace io {

    
    template<>
    void PrettyMatrix<double>::read(jsx::value const& source) {
        data_.resize(0); I_ = source.size(); J_ = I_ ? source(0).size() : 0;
        
        for(auto const& row : source.array()) {
            if(row.size() != J_) throw std::runtime_error("io::prettyrmat: wrong format");
            
            for(auto const& elem : row.array()) data_.push_back(elem.real64());
        }
    };

    template<>
    void PrettyMatrix<double>::write(jsx::value& dest) const {
        dest = jsx::array_t(I_);
        for(int i = 0; i < I_; ++i)
        dest(i) = jsx::array_t(data_.begin() + i*J_, data_.begin() + (i + 1)*J_);
    };


    template<>
    void PrettyMatrix<std::complex<double>>::read(jsx::value const& source) {
        if(!(source.is("real") && source.is("imag") && source.size() == 2))
            throw std::runtime_error("io::prettycmat: wrong format");
        
        data_.resize(0); I_ = source("real").size(); J_ = I_ ? source("real")(0).size() : 0;
        
        if(I_ != source("imag").size())
            throw std::runtime_error("io::prettycmat: wrong format");
        
        for(int i = 0; i < I_; ++i) {
            if(source("real")(i).size() != J_) throw std::runtime_error("io::PrettyMatrix: wrong format");
            if(source("imag")(i).size() != J_) throw std::runtime_error("io::PrettyMatrix: wrong format");
            
            for(int j = 0; j < J_; ++j)
            data_.push_back({source("real")(i)(j).real64(), source("imag")(i)(j).real64()});
        }
    };

    template<>
    void PrettyMatrix<std::complex<double>>::write(jsx::value& dest) const {
        jsx::value real = jsx::array_t(I_, jsx::array_t(J_));
        for(int i = 0; i < I_; ++i)
        for(int j = 0; j < J_; ++j)
        real(i)(j) = operator()(i, j).real();
        
        jsx::value imag = jsx::array_t(I_, jsx::array_t(J_));
        for(int i = 0; i < I_; ++i)
        for(int j = 0; j < J_; ++j)
        imag(i)(j) = operator()(i, j).imag();
        
        dest = jsx::object_t{{"real", real}, {"imag", imag}};
    };
    
};
