#include "Functions.h"

namespace evalsim {
    
    namespace partition {
        
        namespace func {
            
            
            iOmega::iOmega(double beta) : beta_(beta) {};
            std::complex<double> iOmega::operator()(int n) const {
                return {.0, (2*n + 1)*M_PI/beta_};
            };
        
            
            template<typename Value>
            io::Matrix<Value> get_matrix(std::map<std::string, Value> const& entries, jsx::value const& jMatrix)
            {
                io::Matrix<Value> matrix(jMatrix.size(), jMatrix.size());
                
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                        auto const entry = jMatrix(i)(j).string();
                        
                        if(entry != "") matrix(i, j) = entries.at(entry);
                    }
                
                //make symmetric if real
                
                return matrix;
            }
            
            template io::Matrix<double> get_matrix(std::map<std::string, double> const& entries, jsx::value const& jMatrix);
            template io::Matrix<ut::complex> get_matrix(std::map<std::string, ut::complex> const& entries, jsx::value const& jMatrix);
            
            
            template<typename Value>
            std::vector<io::Matrix<Value>> get_function_matrix(std::map<std::string, io::Vector<Value>> const& functions, jsx::value const& jMatrix)
            {
                std::size_t size = std::numeric_limits<std::size_t>::max();
                
                for(auto const& function : functions)
                    size = std::min(size, function.second.size());
                
                std::vector<io::Matrix<Value>> functionMatrix;
                for(std::size_t n = 0; n < size; ++n) {
                    std::map<std::string, Value> entries;
                    
                    for(auto& function : functions)
                        entries[function.first] = function.second[n];
                    
                    functionMatrix.push_back(get_matrix(entries, jMatrix));
                }
                
                return functionMatrix;
            }
            
            template std::vector<io::Matrix<double>> get_function_matrix(std::map<std::string, io::Vector<double>> const& functions, jsx::value const& jMatrix);
            template std::vector<io::Matrix<ut::complex>> get_function_matrix(std::map<std::string, io::Vector<ut::complex>> const& functions, jsx::value const& jMatrix);
            
            
            template<typename Value>
            std::map<std::string, Value> get_entries(io::Matrix<Value> const& matrix, jsx::value const& jMatrix)
            {
                std::map<std::string, std::pair<int, Value>> temp;
                
                for(std::size_t i = 0; i < jMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                        auto const entry = jMatrix(i)(j).string();
                        
                        if(entry != "") {
                            if(!temp.count(entry))
                                temp[entry] = std::pair<int, Value>(0, .0);
                            
                            temp.at(entry).first  += 1;
                            temp.at(entry).second += matrix(i, j);
                        }
                    }
                
                std::map<std::string, Value> entries;
                for(auto const& entry : temp)
                    entries[entry.first] = entry.second.second/Value(entry.second.first);
                
                return entries;
            }
            
            template std::map<std::string, double> get_entries<double>(io::Matrix<double> const& matrix, jsx::value const& jMatrix);
            template std::map<std::string, ut::complex> get_entries<ut::complex>(io::Matrix<ut::complex> const& matrix, jsx::value const& jMatrix);
            
            
            template<typename Value>
            std::map<std::string, io::Vector<Value>> get_function_entries(std::vector<io::Matrix<Value>> const& functionMatrix, jsx::value const& jMatrix)
            {
                std::map<std::string, io::Vector<Value>> functionEntries;
                
                for(auto const& matrix : functionMatrix) {
                    std::map<std::string, Value> entries = get_entries(matrix, jMatrix);
                    for(auto const& entry : entries)
                        functionEntries[entry.first].push_back(entry.second);
                }
                
                return functionEntries;
            }
            
            template std::map<std::string, io::Vector<double>> get_function_entries(std::vector<io::Matrix<double>> const& functionMatrix, jsx::value const& jMatrix);
            template std::map<std::string, io::Vector<ut::complex>> get_function_entries(std::vector<io::Matrix<ut::complex>> const& functionMatrix, jsx::value const& jMatrix);

            
            template<typename Value>
            std::vector<io::Matrix<Value>> get_hybridisation_moment(std::map<std::string, io::cvec> const& functions, jsx::value const& jParams, jsx::value const& jHybMatrix)
            {
                std::map<std::string, io::Vector<Value>> moments;
                io::Vector<Value> tmp(2);
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                        if(jHybMatrix(i)(j).string() != "") {
                            if(jHybMatrix(j)(i).string() == "")
                                throw std::runtime_error("Hyb: invalid hybridisation matrix.");
                            
                            std::string const entry = jHybMatrix(i)(j).string();
                            std::string const entryTransp = jHybMatrix(j)(i).string();
                            
                            if(!moments.count(entry)){
                                auto const fit = bath::Fit<Value>(jParams("beta").real64(), functions.at(entry), functions.at(entryTransp));
                                tmp[0] = fit.moment();
                                tmp[1] = fit.eps();
                                moments[entry] = tmp;
                            }
                        }
                
                return  { get_function_matrix(moments, jHybMatrix) };
            }
            
            template std::vector<io::Matrix<double>> get_hybridisation_moment<double>(std::map<std::string, io::cvec> const& functions, jsx::value const& jParams, jsx::value const& jHybMatrix);
            template std::vector<io::Matrix<ut::complex>> get_hybridisation_moment<ut::complex>(std::map<std::string, io::cvec> const& functions, jsx::value const& jParams, jsx::value const& jHybMatrix);
            
            
            template<typename Value>
            std::tuple<std::vector<io::cmat>, std::vector<io::Matrix<Value>>> get_hybridisation(jsx::value const& jParams)
            {
                jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                jsx::value jFunctions = mpi::read(jParams("hybridisation")("functions").string());
                
                std::map<std::string, io::cvec> functions;
                for(auto& function : jFunctions.object())
                    functions[function.first] = jsx::at<io::cvec>(function.second);
                
                return std::make_tuple(get_function_matrix(functions, jHybMatrix),
                                       get_hybridisation_moment<Value>(functions, jParams, jHybMatrix));  // checks if hybridisation is compatible with Value
            }
            
            template std::tuple<std::vector<io::cmat>, std::vector<io::Matrix<double>>> get_hybridisation<double>(jsx::value const& jParams);
            template std::tuple<std::vector<io::cmat>, std::vector<io::Matrix<ut::complex>>> get_hybridisation<ut::complex>(jsx::value const& jParams);
            
            
            template<typename Value>
            std::vector<io::cmat> get_self_dyson(jsx::value const& jParams, std::vector<io::cmat> const& green, std::vector<io::cmat> const& hyb)
            {
                double const mu = jParams("mu").real64();
                iOmega const iomega(jParams("beta").real64());
                auto const& oneBody = jsx::at<io::Matrix<Value>>(jParams("hloc")("one body"));
                jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                
                auto const size = std::min(green.size(), hyb.size());
                
                std::vector<io::cmat> selfenergy(green.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
                
                for(std::size_t n = 0; n < size; ++n) {
                    io::cmat green_inv = linalg::inv(green[n]);
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                            selfenergy[n](i, j) = (i == j ? iomega(n) + mu : .0) - oneBody(i, j) - hyb[n](i, j) - green_inv(i, j);
                        
                }
                
                return selfenergy;
            }
            
            template std::vector<io::cmat> get_self_dyson<double>(jsx::value const& jParams, std::vector<io::cmat> const& green, std::vector<io::cmat> const& hyb);
            template std::vector<io::cmat> get_self_dyson<ut::complex>(jsx::value const& jParams, std::vector<io::cmat> const& green, std::vector<io::cmat> const& hyb);
            
            
            template<typename Value>
            void add_self_tail(jsx::value const& jParams, std::vector<io::cmat>& function, std::vector<io::Matrix<Value>> const& moments, std::size_t hybSize)
            {
                iOmega const iomega(jParams("beta").real64());
                jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                
                io::Matrix<Value> alpha(jHybMatrix.size(), jHybMatrix.size());
                io::Matrix<Value> gamma(jHybMatrix.size(), jHybMatrix.size());
                
                std::size_t const nFit  = std::max(function.size()/8, static_cast<std::size_t>(1));
                
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                        if(jHybMatrix(i)(j).string() != "") {
                            std::complex<double> iomegaAvg = .0;
                            std::complex<double> selfPosAvg = .0;
                            std::complex<double> selfNegAvg = .0;
                            
                            for(std::size_t n = function.size() - nFit;  n < function.size(); ++n) {
                                iomegaAvg += iomega(n);
                                selfPosAvg += function[n](i, j);
                                selfNegAvg += std::conj(function[n](j, i));
                            }
                            
                            iomegaAvg /= static_cast<double>(nFit);
                            selfPosAvg /= static_cast<double>(nFit);
                            selfNegAvg /= static_cast<double>(nFit);
                            
                            alpha(i, j) = ut::to_val<Value>(-moments[1](i, j)/2.*(1./(selfPosAvg - moments[0](i, j)) + 1./(selfNegAvg - moments[0](i, j))));
                            gamma(i, j) = ut::to_val<Value>(iomegaAvg*moments[1](i, j)/2.*(-1./(selfPosAvg - moments[0](i, j)) + 1./(selfNegAvg - moments[0](i, j))) + iomegaAvg*iomegaAvg);
                        }
                
                
                for(std::size_t n = function.size(); n < hybSize; ++n) {
                    io::cmat temp(jHybMatrix.size(), jHybMatrix.size());
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                            if(jHybMatrix(i)(j).string() != "")
                                temp(i, j) = moments[0](i, j) + moments[1](i, j)/(iomega(n) - alpha(i, j) - gamma(i, j)/iomega(n));
                    
                    function.push_back(temp);
                }
            }
            
            template void add_self_tail(jsx::value const& jParams, std::vector<io::cmat>& function, std::vector<io::Matrix<double>> const& moments, std::size_t hybSize);
            template void add_self_tail(jsx::value const& jParams, std::vector<io::cmat>& function, std::vector<io::Matrix<ut::complex>> const& moments, std::size_t hybSize);
            
            
            template<typename Value>
            void add_green_tail(jsx::value const& jParams, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& selfenergy, std::vector<io::cmat>& green)
            {
                double const mu = jParams("mu").real64();
                iOmega const iomega(jParams("beta").real64());
                auto const& oneBody = jsx::at<io::Matrix<Value>>(jParams("hloc")("one body"));
                jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                
                for(std::size_t n = green.size(); n < hyb.size(); ++n) {
                    io::cmat green_inv(jHybMatrix.size(), jHybMatrix.size());
                    
                    for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                            green_inv(i, j) = (i == j ? iomega(n) + mu : .0) - oneBody(i, j) - hyb[n](i, j) - selfenergy[n](i, j);
                    
                    green.push_back(linalg::inv(green_inv));
                }
            }
            
            template void add_green_tail<double>(jsx::value const& jParams, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& selfenergy, std::vector<io::cmat>& green);
            template void add_green_tail<ut::complex>(jsx::value const& jParams, std::vector<io::cmat> const& hyb, std::vector<io::cmat> const& selfenergy, std::vector<io::cmat>& green);
            
            
            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix, std::vector<io::Matrix<Value>> const& momentsMatrix)
            {
                jsx::value jFunction;
                
                jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                
                std::map<std::string, io::cvec> functions = get_function_entries(functionsMatrix, jHybMatrix);
                for(auto& function : functions)
                    jFunction[function.first]["function"] = std::move(function.second);
                
                std::map<std::string, io::Vector<Value>> moments = get_function_entries(momentsMatrix, jHybMatrix);
                
                for(auto& moment : moments) jFunction[moment.first]["moments"] = moment.second;

                return jFunction;
            }
        
            template jsx::value write_functions(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix, std::vector<io::Matrix<double>> const& momentsMatrix);
            template jsx::value write_functions(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix, std::vector<io::Matrix<ut::complex>> const& momentsMatrix);
            
            
            template<typename Value>
            jsx::value write_functions(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix)
            {
                jsx::value jFunction;
                
                jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
                
                std::map<std::string, io::cvec> functions = get_function_entries(functionsMatrix, jHybMatrix);
                for(auto& function : functions)
                    jFunction[function.first]["function"] = std::move(function.second);
                
                return jFunction;
            }
            
            template jsx::value write_functions<double>(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix);
            template jsx::value write_functions<ut::complex>(jsx::value const& jParams, std::vector<io::cmat> const& functionsMatrix);
         
            
            
            template<typename Value>
            std::vector<io::cmat> get_aux_green(jsx::value const& jParams, std::vector<io::cmat> const& selfenergy, std::vector<io::Matrix<Value>> const& selfMoments)
            {
                
                iOmega const iomega(jParams("beta").real64());
                jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
                
                auto const size = selfenergy.size();
                
                std::size_t const norb = jHybMatrix.size();
                
                std::vector<io::cmat> aux(selfenergy.size(), io::cmat(norb, norb));
                
                for(std::size_t n = 0; n < size; ++n){
                    for(std::size_t i = 0; i < norb; ++i)
                        aux[n](i, i) = 1./(iomega(n) + selfMoments[0](i,i)  - selfenergy[n](i, i));
                        
                    for(std::size_t i = 0; i < norb; ++i)
                        for(std::size_t j = 0; j < norb; ++i)
                            if (i!=j){
                                aux[n](i, j) = -2./(iomega(n) + selfMoments[0](i,j)  - selfenergy[n](i, j));
                                aux[n](i, j) += aux[n](i, i) + aux[n](j, j);
                                aux[n](i, j) *= 0.5;
                            }
                    
                }
                
                return aux;
            }
            
            template std::vector<io::cmat> get_aux_green<double>(jsx::value const& jParams, std::vector<io::cmat> const& selfenergy, std::vector<io::Matrix<double>> const& selfMoments);
            template std::vector<io::cmat> get_aux_green<ut::complex>(jsx::value const& jParams, std::vector<io::cmat> const& selfenergy, std::vector<io::Matrix<ut::complex>> const& selfMoments);
            
            
            jsx::value fourier_transform(jsx::value const& jParams, jsx::value const& jGm, int const nf, int const nhf, int const ntau){
                
                auto const beta = jParams("beta").real64();
                
                auto const taus = linalg::linspace<double>(ntau, 0, beta);
                
                iOmega const iomega(beta);
                std::vector<double> omegas(nf);
                for (int i=0; i < omegas.size(); i++)
                    omegas[i] = iomega(i).imag();
                
                jsx::value jGtau;
                for (auto& function : jGm.object()){
                    auto f = function.second;
                    auto Gm = jsx::at<io::cvec>(f("function"));
                    jGtau[function.first] = io::cvec(ft::inverseFourier( Gm, omegas, taus, beta, nhf));
                }
                
                return jGtau;
                
            }
            
            
            
        }
        
    }
    
}









