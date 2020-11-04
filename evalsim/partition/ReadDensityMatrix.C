#include "ReadDensityMatrix.h"


namespace evalsim {
    
    namespace partition {
        
        namespace meas {

            template<typename Value> // one coud get rid of this fuck when properly implementing meas::Matrix<Value> ....
            jsx::value read_matrix(io::Vector<Value> const& source, int const dim) {
                if(dim*dim != source.size())
                    throw std::runtime_error("read_matrix: missmatch in size");
                
                io::Matrix<Value> dest(dim, dim);
                for(std::size_t i = 0; i < dim; ++i)
                    for(std::size_t j = 0; j < dim; ++j)
                        dest(i, j) = source.at(i*dim + j);
                
                return dest;
                
            };
            
            template jsx::value read_matrix<double>(io::Vector<double> const& source, int const dim) ;
            template jsx::value read_matrix<ut::complex>(io::Vector<ut::complex> const& source, int const dim) ;
            
            template<typename Value>
            jsx::value read_density_matrix(jsx::value const& jParams, jsx::value const& jDensityMatrixMeas) {
                jsx::value const& jEigenValues = jParams("hloc")("eigen values");
                
                jsx::value jDensityMatrix = jsx::array_t(jEigenValues.size());
                
                for(std::size_t sec = 0; sec < jEigenValues.size(); ++sec) {
                    jDensityMatrix(sec)["target"] = jsx::int64_t(sec);
                    jDensityMatrix(sec)["matrix"] = read_matrix(jsx::at<io::Vector<Value>>(jDensityMatrixMeas(sec)), jsx::at<io::rvec>(jEigenValues(sec)).size());
                }
                
                return jDensityMatrix;
            };
            
            template jsx::value read_density_matrix<double>(jsx::value const& jParams, jsx::value const& jDensityMatrixMeas);
            template jsx::value read_density_matrix<ut::complex>(jsx::value const& jParams, jsx::value const& jDensityMatrixMeas);
            
            
        }
        
    }
    
}

