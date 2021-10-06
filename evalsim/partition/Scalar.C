#include "Scalar.h"

namespace evalsim {
    
    namespace partition {
        
        template<typename Value>
        jsx::value get_scalar(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements)
        {
            jsx::value jScalar;

            jScalar["k"] = jMeasurements("scalar")("k");
            
            mpi::cout << "Calculating quantum number observables ... " << std::flush;

            auto const& sectorProb = jsx::at<io::rvec>(jMeasurements("sector prob"));
            
            for(auto const& jqn : jPartition("quantum numbers").object()) {
                auto Qn = ga::construct_sector_qn(jParams("hloc"), jqn.second);
                
                double val = .0, valval = .0;
                for(int s = 0; s < jParams("hloc")("eigen values").size(); ++s) {
                    val += sectorProb.at(s)*Qn.at(s);
                    valval += sectorProb.at(s)*Qn.at(s)*Qn.at(s);
                }
                
                jScalar[jqn.first] = io::rvec{{val}};
                jScalar[jqn.first + jqn.first] = io::rvec{{valval}};
            }
            
            mpi::cout << "Ok" << std::endl;
            
            
            jsx::value jDensityMatrix = meas::read_density_matrix<Value>(jParams, jMeasurements("density matrix"));
            

            mpi::cout << "Calculating observables ... " << std::flush;
            
            for(auto const& jObs : jPartition("observables").object())
                jScalar[jObs.first] = io::rvec{{ut::real(linalg::trace<Value>(jDensityMatrix, jObs.second))}};
            
            mpi::cout << "Ok" << std::endl;
            
            
            mpi::cout << "Calculating energy ... " << std::flush;
             
            jsx::value jHamiltonianEff = get_effective_hamiltonian<Value>(jParams);
            
            double energy = ut::real(linalg::trace<Value>(jDensityMatrix, jHamiltonianEff));
            if(jParams.is("dyn")) energy += jsx::at<io::rvec>(jMeasurements("dyn energy"))[0];
            jScalar["energy"] = io::rvec{{energy}};
            
            mpi::cout << "Ok" << std::endl;
            
            
            mpi::cout << "Calculating (log of) the Atomic Partition Function ... " << std::flush;
            
            double beta = jParams("beta").real64();
            io::rvec energies;
            
            for(std::size_t s = 0; s < jHamiltonianEff.size(); ++s){
                if(!jHamiltonianEff(s)("target").is<jsx::null_t>()) {
                    std::size_t const m = jHamiltonianEff(s)("target").int64();
                    
                    if(!jHamiltonianEff(m)("target").is<jsx::null_t>()) {
                        std::size_t const t = jHamiltonianEff(m)("target").int64();
                        
                        auto mat = jsx::at<io::rmat>(jHamiltonianEff(s)("matrix"));
                        for (std::size_t i = 0; i < mat.I(); ++i){
                            energies.push_back( - beta * mat(i,i));
                        }
                    }
                }
            }
            
            auto it = std::max_element(std::begin(energies), std::end(energies));
            auto const maxe = *it;
            auto sum = 0;
            for (auto const& x : energies)
                sum += std::exp(x - maxe);
            
            jScalar["ln(atomic Z)"] = io::rvec{{maxe + std::log(sum)}};
                
                
            mpi::cout << "Ok" << std::endl;
            
            return jScalar;
        }
        
        template jsx::value get_scalar<double>(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);
        template jsx::value get_scalar<ut::complex>(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements);

        
    }
}









