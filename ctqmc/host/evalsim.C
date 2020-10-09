

#include "../../evalsim/Evalsim.h"

namespace evalsim {
    
    void add_dynamics(jsx::value const& jParams, jsx::value& jMeasurements, std::string const worm, std::string const meas)
    {
        if(jParams.is("dyn") && jParams.is(worm + meas)) {
            auto& jDynamic = jMeasurements(worm)("dynamic");
            auto& jStatic = jMeasurements(worm + meas)("static");

            std::set<std::string> entries;
        
            for(auto entry : jDynamic.object()) entries.insert(entry.first);
            for(auto entry : jStatic.object())  entries.insert(entry.first);
        
            for(auto entry : entries) {
                if(jDynamic.is(entry) && !jStatic.is(entry))
                    jStatic[entry] = jDynamic(entry);
                else if (jDynamic.is(entry) && jStatic.is(entry)) {
                    if(!add_dynamic<double>(jStatic(entry), jDynamic(entry)) && !add_dynamic<ut::complex>(jStatic(entry), jDynamic(entry)))
                        throw std::runtime_error("evalsim::add_dynamics: missmatch for " + worm + meas);
                }
            }
        }
    }
    
    
    namespace partition {
        
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
        
        double serial_number(int sector, int i){
            return std::stod(std::to_string(sector) + "." + std::to_string(i));
        }
        
        jsx::value get_qn_susc(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jMeasurements, jsx::value const& jScalar) {
            double const beta = jParams("beta").real64(); jsx::value const& jHybMatrix = jParams("hybridisation")("matrix");
            
            jsx::value jSusc;
            
            for(auto jqn : jPartition("quantum numbers").object()) {
                
                mpi::cout << "Reading " << jqn.first << " susceptibility ... " << std::flush;
                
                auto const qn = jsx::at<io::rvec>(jqn.second);
                
                double moment = 0;
                io::rvec function(jsx::at<io::rvec>(jMeasurements("susceptibility flavor")("0_0")).size(), .0);
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                        double const fact = qn.at(i)*qn.at(j);
                        auto const meas = jsx::at<io::rvec>(jMeasurements("susceptibility flavor")(std::to_string(i) + "_" + std::to_string(j)));
                        
                        function[0] += fact*meas[0]/(2.*beta);
                        
                        for(std::size_t n = 1; n < function.size(); ++n) function[n] += fact*beta*meas[n]/(4*M_PI*M_PI*n*n);
                        
                        if(i == j) moment -= fact*(jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;
                    }
                function[0] += beta*jsx::at<io::rvec>(jScalar(jqn.first + jqn.first)).at(0);
                function[0] -= beta*jsx::at<io::rvec>(jScalar(jqn.first)).at(0)*jsx::at<io::rvec>(jScalar(jqn.first)).at(0);
                
                mpi::cout << "Ok" << std::endl;
                
                
                mpi::cout << "Adding " << jqn.first << " susceptibility high frequency tail ... " << std::flush;
                
                std::size_t const nFit = std::max(static_cast<int>(function.size()/10.), 1);
                
                double A = .0, B = .0;
                for(std::size_t n = function.size() - nFit; n < function.size(); ++n) {
                    A += moment*function[n] + (4*M_PI*M_PI*n*n/(beta*beta))*function[n]*function[n]; B += function[n]*function[n];
                }
                double const alpha = -A/B;
                
                std::size_t const nTail = std::max(static_cast<int>(beta*jPartition("susceptibility tail").real64()/(2*M_PI)), 1);
                for(std::size_t n = function.size(); n < nTail; ++n)
                    function.push_back(-moment/((4*M_PI*M_PI*n*n/(beta*beta)) + alpha));
                
                jSusc[jqn.first]["function"] = function;
                jSusc[jqn.first]["moment"] = io::rvec{{ moment }};
                
                mpi::cout << "Ok" << std::endl;
            }
            
            
            return jSusc;
        }
        
        namespace meas{
            
            namespace impl{
                
                std::vector<io::cmat> symmetrize_functions(std::vector<io::cmat> const& functions, std::vector<io::cmat> const& functionsConj, jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jHybMatrix, int hybSize)
                {
                    std::vector<io::cmat> symmetrized(functions.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
                    
                    for(std::size_t n = 0; n < symmetrized.size(); ++n)
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                            for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                                symmetrized[n](i, j) = (functions[n](i, j) + std::conj(functionsConj[n](j, i)))/2.;

                    return symmetrized;
                };
                
                
            }
            
        }
        
    }
    
    namespace worm {
        
        std::string const Kernels::name = "kernels";
        
        std::vector<io::ctens> rearrange_susc_ph(std::vector<io::ctens> const& susc_ph){
            std::vector<io::ctens> r(susc_ph.size(),io::ctens(susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I()));
            for (auto const& ijkl : susc_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //j and i are swapped
                //k and l are swapped
                std::string const entry = std::to_string(2*j)+"_"+std::to_string(2*i+1)+"_"+std::to_string(2*l)+"_"+std::to_string(2*k+1);
                
                for (int n=0; n<susc_ph.size(); n++)
                    r[n].emplace(j,i,l,k, entry, -susc_ph[n].at(i,j,k,l));
                
            }

            return r;
        }

        std::vector<io::ctens> rearrange_susc_tph(std::vector<io::ctens> const& susc_ph){
            std::vector<io::ctens> r(susc_ph.size(),io::ctens(susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I(),susc_ph[0].I()));
            
            for (auto const& ijkl : susc_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //j and l are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*l+1)+"_"+std::to_string(2*k)+"_"+std::to_string(2*j+1);
                
                for (int n=0; n<susc_ph.size(); n++)
                    r[n].emplace(i,l,k,j, entry, -susc_ph[n].at(i,j,k,l));
                
            }

            return r;
            
        }

        std::vector<io::ctens> rearrange_susc_pp(std::vector<io::ctens> const& susc_pp){
            std::vector<io::ctens> r(susc_pp.size(),io::ctens(susc_pp[0].I(),susc_pp[0].I(),susc_pp[0].I(),susc_pp[0].I()));

            for (auto const& ijkl : susc_pp[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //j and k are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*k+1)+"_"+std::to_string(2*j)+"_"+std::to_string(2*l+1);
                
                for (int n=0; n<susc_pp.size(); n++)
                    r[n].emplace(i,k,j,l, entry, -susc_pp[n].at(i,j,k,l));
                
            }

            return r;

        }


        std::vector<io::ctens> rearrange_hedin_ph(std::vector<io::ctens> const& hedin_ph){
            std::vector<io::ctens> r(hedin_ph.size(),io::ctens(hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I()));
            
            for (auto const& ijkl : hedin_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //k and l are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*j+1)+"_"+std::to_string(2*l)+"_"+std::to_string(2*k+1);

                for (int n=0; n<hedin_ph.size(); n++){
                    r[n].emplace(i,j,l,k, entry, hedin_ph[n].at(i,j,k,l));
                }
                
            }

            return r;

        }

        std::vector<io::ctens> rearrange_hedin_tph(std::vector<io::ctens> const& hedin_ph){
            std::vector<io::ctens> r(hedin_ph.size(),io::ctens(hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I(),hedin_ph[0].I()));

            for (auto const& ijkl : hedin_ph[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //l and j are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*l+1)+"_"+std::to_string(2*k)+"_"+std::to_string(2*j+1);
                
                for (int n=0; n<hedin_ph.size(); n++)
                    r[n].emplace(i,l,k,j, entry, hedin_ph[n].at(i,j,k,l));
                
            }

            return r;
        }

        std::vector<io::ctens> rearrange_hedin_pp(std::vector<io::ctens> const& hedin_pp){
            std::vector<io::ctens> r(hedin_pp.size(),io::ctens(hedin_pp[0].I(),hedin_pp[0].I(),hedin_pp[0].I(),hedin_pp[0].I()));

            for (auto const& ijkl : hedin_pp[0].ijkl()){

                auto const i=ijkl[1];
                auto const j=ijkl[2];
                auto const k=ijkl[3];
                auto const l=ijkl[4];
                
                //k and j are swapped
                std::string const entry = std::to_string(2*i)+"_"+std::to_string(2*k+1)+"_"+std::to_string(2*j)+"_"+std::to_string(2*l+1);
                                   
                for (int n=0; n<hedin_pp.size(); n++)
                    r[n].emplace(i,k,j,l, entry, -hedin_pp[n].at(i,j,k,l));
                
                
            }

            return r;
        }
        
        namespace meas {
            
            namespace impl {
                
                io::cvec read_matsubara_function(jsx::value const& jParams, jsx::value const& jPartition, jsx::value const& jFunction, bool scale){
                    double const beta = jParams("beta").real64();
                    
                    io::cvec function = jsx::at<io::cvec>(jFunction);
                    
                    if (scale) for (auto& x : function) x*=-1./beta;
                    
                    return function;
                }
                
                io::cvec collapse_antisymmetric_function(io::cvec const& function){
                    
                    io::cvec as_function(function.size()/2);
                    
                    std::size_t start = function.size()/2;
                    
                    auto it_forward = function.begin()+start;
                    auto it_backward = function.begin()+start-1;
                    
                    for (auto& x : as_function){
                        x = 0.5 * ( *it_forward + std::conj(*it_backward) );
                        --it_backward; ++it_forward;
                    }
                    
                    return as_function;
                }
                
            }
            
        }
        
        namespace func {
            
            std::vector<io::cmat> green_function_on_full_axis(std::vector<io::cmat>const& green){
                std::vector<io::cmat> r(2*green.size(),io::cmat(green[0].I(),green[0].J()));
                
                auto it_forward = r.begin() + green.size();
                auto it_backward = r.begin() + green.size()-1;
                for(auto const& x : green){
                    *it_forward++ = x;
                    *it_backward-- = x.conj();
                }
                    
                return r;
            }
            
        }
        
    }
}
