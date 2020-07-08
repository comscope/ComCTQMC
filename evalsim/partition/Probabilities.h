#ifndef EVALSIM_PARTITION_PROBABILITIES_H
#define EVALSIM_PARTITION_PROBABILITIES_H


#include "ReadDensityMatrix.h"
#include "ReadHamiltonian.h"


#include "../../include/linalg/Operators.h"
#include "../../include/JsonX.h"
#include "../../include/io/Vector.h"
#include "../../include/io/Matrix.h"
#include "../../include/options/Options.h"
#include "../../include/atomic/Generate.h"

namespace evalsim {
    
    namespace partition {
        
        
        struct VecLess {
            VecLess(std::size_t size, std::size_t start, std::size_t end) : size_(size), start_(start), end_(end) {};
            bool operator()(std::vector<double> const& lhs, std::vector<double> const& rhs) const {
                if(lhs.size() != size_ || rhs.size() != size_)
                    throw std::runtime_error("VecCompare");
                
                for(std::size_t i = start_; i < end_; ++i) {
                    if(lhs[i] > rhs[i])
                        return true;
                    else if(lhs[i] < rhs[i])
                        return false;
                }
                
                return false;
            }
        private:
            std::size_t const size_, start_, end_;
        };
        
        
        double truncate(double val, int prec) {
            if(std::abs(val) < 1.e-8) return .0;
            
            std::stringstream temp;
            temp << std::setprecision(prec) << val;
            temp >> val; return val;
        }
    
        ut::complex truncate(ut::complex val, int prec) {
            if(std::abs(val.real()) < 1.e-8) val.real(0.);
            if(std::abs(val.imag()) < 1.e-8) val.imag(0.);
            
            std::stringstream temp;
            temp << std::setprecision(prec) << val;
            temp >> val; return val;
        }

        template<typename Value>
        std::vector<std::string> get_surviving_qn(jsx::value const& jPartition){
            
            std::vector<std::string> surviving;
            for(int qn = 0; qn < jPartition("probabilities").size(); ++qn) {
                std::string const name = jPartition("probabilities")(qn).string();
                if(jPartition("quantum numbers").is(name))
                    surviving.push_back(name);
            }
    
            return surviving;
        }
        
        template<typename Value>
        jsx::array_t get_eigenstate_probabilities(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements){
            
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            
            std::vector<std::vector<double>> data;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector){
                
                auto const N = jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size();
                std::vector<double> temp(N);
                
                for(int i = 0; i < N; ++i) {
                    temp[i] = std::abs(jsx::at<io::Matrix<Value>>(jDensityMatrix(sector)("matrix"))(i, i));     // take abs(real) value ?
                   }
                
                data.push_back(temp);
            }
            
            jsx::array_t result;
            for(auto& entry : data) {
                result.push_back(jsx::array_t(entry.begin(), entry.end()));
            }
            
            return result;
        }
    
        template<typename Value>
        jsx::array_t get_eigenstate_energies(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements){
            
            jsx::value jHamiltonianEff = get_effective_hamiltonian<Value>(jParams);
            
            std::vector<std::vector<double>> data;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector){
                
                auto const N = jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size();
                std::vector<double> temp(N);
                
                for(int i = 0; i < N; ++i) {
                    temp[i] = truncate(ut::real(jsx::at<io::Matrix<Value>>(jHamiltonianEff(sector)("matrix"))(i, i)), 8);
                   }
                
                data.push_back(temp);
            }
            
            jsx::array_t result;
            for(auto& entry : data) {
                result.push_back(jsx::array_t(entry.begin(), entry.end()));
            }
            
            return result;
        }
    
        template<typename Value>
        jsx::array_t get_eigenstate_qn(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements){
            
            auto const surviving = get_surviving_qn<Value>(jPartition);
            
            std::vector<std::vector<double>> data;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i){
                    std::vector<double> temp(surviving.size());
                    data.push_back(temp);
                }
            
            for(int qn = 0; qn < surviving.size(); ++qn) {
                auto const name = surviving[qn];
                
                if(jPartition("quantum numbers").is(name)) {
                    auto Qn = ga::construct_sector_qn(jParams("hloc"), jPartition("quantum numbers")(name));
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(Qn.at(sector), 8);
                    
                }
            }
            
            jsx::array_t result;
            for(auto& entry : data) {
                result.push_back(jsx::array_t(entry.begin(), entry.end()));
            }
            
            return result;
        }
        
        
        template<typename Value>
        jsx::value get_eigenstates(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements)
        {
            jsx::value jOverlap;
            
            mpi::cout << "describing eigenstates ... " << std::flush;
            
            jsx::value jTransformation = jParams("hloc")("transformation"); // This looks like the eigenvectors of hloc (sector)("matrix")(i,j)
            jsx::value jOccupationStates = ga::construct_occupation_states<Value>(jParams("hloc"));
            
            auto const qn = get_eigenstate_qn<Value>(jParams, jPartition, jMeasurements);
            auto const prob = get_eigenstate_probabilities<Value>(jParams, jPartition, jMeasurements);
            auto const energies = get_eigenstate_energies<Value>(jParams, jPartition, jMeasurements);
            auto const surviving = get_surviving_qn<Value>(jPartition);
            jOverlap["surviving"] = jsx::array_t(surviving.begin(),surviving.end());
            
            int const n_orb = jParams("hybridisation")("matrix").array().size();
            jsx::value data; int index = 0;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector){
                auto const eigenvectors = jsx::at<io::Matrix<Value>>(jTransformation(sector)("matrix"));
                int const N = jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size();
                
                io::Matrix<int> occupation_basis(n_orb,N);
                io::Matrix<Value> occupation_representations(n_orb,N);
                io::Matrix<Value> results(n_orb,N);
                auto occupation_representations_ptr = occupation_representations.data();
                auto occupation_basis_ptr = occupation_basis.data();
                
                //std::vector<std::string> basis(N,"");
                for(int i = 0; i < N ; ++i) {
                    auto const occupation_representation = jsx::at<io::rvec>(jOccupationStates(index++));
                    for (auto const& x : occupation_representation){
                        *occupation_representations_ptr++ = x;
                        *occupation_basis_ptr++ = x;
                        //basis[i]+=std::to_string(static_cast<int>(x));
                    }
                }
                
                /* testing
                if (sector==100){
                    std::cout << "\n";
                    for (auto& x : basis) std::cout << x << "\n";
                    for (int i = 0; i < eigenvectors.I(); i++){
                        for (int j = 0; j < eigenvectors.J(); j++)
                            std::cout << ut::real(eigenvectors(i,j)) << " ";
                        std::cout << "\n";
                    }
                }
                */
                
                linalg::mult<Value>('n', 't', 1., occupation_representations, eigenvectors, .0, results);
                
                for (int i = 0; i < results.I(); i++)
                    for (int j = 0; j < results.J(); j++)
                        results(i,j) = truncate(results(i,j),8);
                
                
                data[std::to_string(sector)]["occupation representation basis"] = std::move(occupation_basis);
                data[std::to_string(sector)]["eigenvectors"] = std::move(eigenvectors);
                data[std::to_string(sector)]["occupation representation"] = std::move(results);
                data[std::to_string(sector)]["quantum numbers"] = qn[index-1]; // should be identical for all index in a sector
                data[std::to_string(sector)]["probabilities"] = prob[sector];
                data[std::to_string(sector)]["energies"] = energies[sector];
                
            }
            
            jOverlap = std::move(data);
            
            mpi::cout << "Ok" << std::endl;
            
            return jOverlap;
        }
        
        
        template<typename Value>
        jsx::value get_probabilities(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements)
        {
            jsx::value jProbabilities;
            
            mpi::cout << "Reading impurity probabilities ... " << std::flush;
            
            jsx::value jHamiltonianEff = get_effective_hamiltonian<Value>(jParams);
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            
            
            //Collect list of probabilities that have survived the criterion
            std::vector<std::string> surviving;
            for(int qn = 0; qn < jPartition("probabilities").size(); ++qn) {
                std::string const name = jPartition("probabilities")(qn).string();
                if(name == "energy" || jPartition("quantum numbers").is(name) || jPartition("observables").is(name))
                    surviving.push_back(name);
            }
            
            jsx::array_t temp;
            for(auto& entry : surviving) {
                temp.push_back(jsx::string_t(entry));
            }
            jProbabilities["surviving"] = std::move(temp);
            
            std::vector<std::vector<double>> data;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i) {
                    std::vector<double> temp(surviving.size() + 1);
                    temp.back() = std::abs(jsx::at<io::Matrix<Value>>(jDensityMatrix(sector)("matrix"))(i, i));     // take abs(real) value ?
                    data.push_back(temp);
                }
            
            for(int qn = 0; qn < surviving.size(); ++qn) {
                std::string const name = surviving[qn];
                
                if(jPartition("quantum numbers").is(name)) {
                    auto Qn = ga::construct_sector_qn(jParams("hloc"), jPartition("quantum numbers")(name));
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(Qn.at(sector), 8);
                    
                } else if(jPartition("observables").is(name)) {
                    jsx::value const& jObservable = jPartition("observables")(name);
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn]= truncate(ut::real(jsx::at<io::Matrix<Value>>(jObservable(sector)("matrix"))(i, i)), 8);
                    
                } else if(name == "energy") {
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(ut::real(jsx::at<io::Matrix<Value>>(jHamiltonianEff(sector)("matrix"))(i, i)), 8);
                    
                }
            }
            
            
            if(surviving.size()) std::sort(data.begin(), data.end(), VecLess(surviving.size() + 1, 0, surviving.size()));
            
            
            temp.clear();
            for(auto& entry : data) {
                temp.push_back(jsx::array_t(entry.begin(), entry.end()));
                temp.back().array().back() = io::rvec{{entry.back()}};
            }
            jProbabilities["quantum numbers"] = std::move(temp);
            
            
            
            {
                jsx::value jTransformation = jParams("hloc")("transformation"), jTemp;
                
                linalg::mult<Value>('n', 'c', 1., jDensityMatrix, jTransformation, .0, jTemp);
                linalg::mult<Value>('n', 'n', 1., jTransformation, jTemp, .0, jDensityMatrix);
            }
            
            jsx::value jOccupationStates = ga::construct_occupation_states<Value>(jParams("hloc"));
            
            data.clear(); int index = 0;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i) {
                    io::rvec temp = jsx::at<io::rvec>(jOccupationStates(index++));
                    temp.push_back(std::abs(jsx::at<io::Matrix<Value>>(jDensityMatrix(sector)("matrix"))(i, i)));       // take abs(real) value ?
                    data.push_back(temp);
                }
            
            jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            
            std::sort(data.begin(), data.end(), VecLess(jHybMatrix.size() + 1, jHybMatrix.size(), jHybMatrix.size() + 1));
            
            temp.clear();
            for(auto& entry : data) {
                temp.push_back(jsx::array_t(entry.begin(), entry.end()));
                temp.back().array().back() = io::rvec{{entry.back()}};
            }
            jProbabilities["occupation numbers"] = std::move(temp);
            
            mpi::cout << "Ok" << std::endl;
            
            return jProbabilities;
        }
   
    }
    
}


#endif









