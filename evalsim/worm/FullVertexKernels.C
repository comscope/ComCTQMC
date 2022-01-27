#include "FullVertexKernels.h"

namespace evalsim {
    
    namespace worm {
        
        std::string const Kernels::name = "kernels";
        
        //We must swap operator orders to get everything into the same order as the full vertex
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
    
        template<typename Value>
        jsx::value evaluateKernels(jsx::value jParams, jsx::value const& jObservables) {
            
            mpi::cout << "Evaluating asymptotic vertex kernels" << std::endl;
            
            jsx::value jWorm = jParams(worm::Kernels::name);
            
            int b_cutoff =  jWorm.is("boson cutoff") ? std::max(1, int(jWorm("boson cutoff").int64())) : 1;
            int f_cutoff = jWorm.is("fermion cutoff") ? std::max(1, int(jWorm("fermion cutoff").int64())) : 50;
            
            
            //double const beta = jParams("beta").real64();
            int const pos_and_neg_freq = (std::is_same<Value,double>::value ? 1 : 2);
            int const nMatGB = pos_and_neg_freq*b_cutoff - pos_and_neg_freq + 1;
            int const nMatGF = 2*f_cutoff;
            
            func::OmegaMap omega_f(nMatGF,false,true);
            func::OmegaMap omega_b(nMatGB,true,!std::is_same<Value,double>::value);
            
            int const nMatGB_kernel = 2*b_cutoff - 1;
            int const nMatGF_kernel = nMatGF;
            
            func::OmegaMap omega_f_kernel(nMatGF_kernel,false,true);
            func::OmegaMap omega_b_kernel(nMatGB_kernel,true,true);
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            //Read in green functions and susc and hedin susceptibilities
            mpi::cout << "Reading in Susceptibilities and Green's functions ... " << std::flush;
            
            auto green = std::vector<io::cmat>(); auto self = std::vector<io::cmat>();
            func::green_and_self_function_on_full_axis<Value>(jParams, jObservables, jHybMatrix, hyb, hybMoments, 2*nMatGB + nMatGF, green, self);
            func::OmegaMap greenOM(green.size(),false,true);
            
            mpi::cout << " ... susc ... " << std::flush;
            
            auto const susc_ph = jObservables.is(cfg::susc_ph::Worm::name()) ?
            rearrange_susc_ph( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::susc_ph::Worm::name())("susceptibility"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size())) : std::vector<io::ctens>();
            auto const susc_pp = jObservables.is(cfg::susc_pp::Worm::name()) ?
                rearrange_susc_pp( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::susc_pp::Worm::name())("susceptibility"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size())) : std::vector<io::ctens>();
            auto const susc_tph = jObservables.is(cfg::susc_ph::Worm::name()) ? rearrange_susc_tph(susc_ph) : std::vector<io::ctens>();
            
            auto const do_sph = susc_ph.size() > 0;
            auto const do_spp = susc_pp.size() > 0;
            auto const do_stph = susc_tph.size() > 0;
            
            mpi::cout << " ... hedin ... " << std::flush;
            
            auto const hedin_ph = jObservables.is(cfg::hedin_ph_imprsum::Worm::name()) ?
                rearrange_hedin_ph( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::hedin_ph_imprsum::Worm::name())("susceptibility"), jParams, jParams(cfg::hedin_ph_imprsum::Worm::name()), jHybMatrix, hyb.size())) : std::vector<io::ctens>();
            auto const hedin_pp = jObservables.is(cfg::hedin_pp_imprsum::Worm::name()) ?
            rearrange_hedin_pp( meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::hedin_pp_imprsum::Worm::name())("susceptibility"), jParams, jParams(cfg::hedin_pp_imprsum::Worm::name()), jHybMatrix, hyb.size())) : std::vector<io::ctens>();
            auto const hedin_tph = jObservables.is(cfg::hedin_ph_imprsum::Worm::name()) ? rearrange_hedin_tph(hedin_ph) : std::vector<io::ctens>();
            
            auto const do_hph = hedin_ph.size() > 0;
            auto const do_hpp = hedin_pp.size() > 0;
            auto const do_htph = hedin_tph.size() > 0;
            
            mpi::cout << " ... vertex ... " << std::flush;
            
            std::vector<io::ctens> vertex = jObservables.is(cfg::vertex_imprsum::Worm::name()) ?
            meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::vertex_imprsum::Worm::name())("susceptibility"), jParams, jParams(cfg::vertex_imprsum::Worm::name()), jHybMatrix, hyb.size()) : std::vector<io::ctens>();
            
            auto const do_v = vertex.size() > 0;
            
            mpi::cout << "OK" << std::endl;
            
            //Construct interaction matrix U_ijkl
            mpi::cout << "Constructing interaction matrix ... " << std::flush;
            
            //This matrix has a factor of 1/2 built in, so we must adjust the resulting kernel equations
            //params::complete_impurity<Value>(jParams);
            imp::Tensor<Value> const U_tmp(jParams("hloc")("two body"),jHybMatrix.size());
            InteractionTensor<Value> const U(U_tmp,jHybMatrix.size());
            
            mpi::cout << "OK" << std::endl;
            
            //Construct Kernel-1 Functions
            mpi::cout << "Calculating Kernel-1 functions ... " << std::flush;
            
            std::vector<io::ctens> kernel_0(1, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_1_ph(nMatGB_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_1_pp(nMatGB_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_1_tph(nMatGB_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            auto const abcds = do_v ? vertex[0].ijkl() : construct_ijkls(jHybMatrix.size(), U);
            
            int const size = nMatGB_kernel;
            int const chunk = (size + mpi::number_of_workers() - 1)/mpi::number_of_workers();
            int const start = chunk*mpi::rank();
            int const end = chunk*(mpi::rank() + 1);
            
            for (int n=start; n<end; n++){
                if (n >= nMatGB_kernel) break;
                
                int const n_susc_test = omega_b.pos(omega_b_kernel(n));
                bool const do_conj = n_susc_test < 0 ? true : false;
                int const n_susc = do_conj ? omega_b.pos(-omega_b_kernel(n)) : n_susc_test;
                
                for(auto const& abcd : abcds){
                    
                    auto const a=abcd[1];
                    auto const b=abcd[2];
                    auto const c=abcd[3];
                    auto const d=abcd[4];
                    
                    std::string const entry = std::to_string(2*a)+"_"+std::to_string(2*b+1)+"_"+std::to_string(2*c)+"_"+std::to_string(2*d+1);
                    
                    if (!n) kernel_0[n].emplace(a,b,c,d,entry,-2.*U(a,c,b,d));
                    
                    kernel_1_ph[n].emplace(a,b,c,d,entry,0.);
                    kernel_1_tph[n].emplace(a,b,c,d,entry,0.);
                    kernel_1_pp[n].emplace(a,b,c,d,entry,0.);
                    
                    if (do_sph)
                        for (auto const& ijkl : susc_ph[0].ijkl()){
                            auto const i=ijkl[1];
                            auto const j=ijkl[2];
                            auto const k=ijkl[3];
                            auto const l=ijkl[4];
                            if (do_conj)
                                kernel_1_ph[n](a,b,c,d) -= 4.*U(a,j,b,i) * std::conj(susc_ph[n_susc].at(i,j,k,l)) * U(l,c,k,d);
                            else
                                kernel_1_ph[n](a,b,c,d) -= 4.*U(a,j,b,i) * susc_ph[n_susc].at(i,j,k,l) * U(l,c,k,d);
                            
                        }
                      
                    if (do_stph)
                        for (auto const& ijkl : susc_tph[0].ijkl()){
                            auto const i=ijkl[1];
                            auto const j=ijkl[2];
                            auto const k=ijkl[3];
                            auto const l=ijkl[4];
                            if (do_conj)
                                kernel_1_tph[n](a,b,c,d) -= 4.*U(a,l,i,d) * std::conj(susc_tph[n_susc].at(i,j,k,l)) * U(j,c,b,k);
                            else
                                kernel_1_tph[n](a,b,c,d) -= 4.*U(a,l,i,d) * susc_tph[n_susc].at(i,j,k,l) * U(j,c,b,k);
                        }
                    
                    if (do_spp)
                        for (auto const& ijkl : susc_pp[0].ijkl()){
                            auto const i=ijkl[1];
                            auto const j=ijkl[2];
                            auto const k=ijkl[3];
                            auto const l=ijkl[4];
                            if (do_conj)
                                kernel_1_pp[n](a,b,c,d) -= U(a,c,k,i) * std::conj(susc_pp[n_susc].at(i,j,k,l)) * U(l,j,b,d);
                            else
                                kernel_1_pp[n](a,b,c,d) -= U(a,c,k,i) * susc_pp[n_susc].at(i,j,k,l) * U(l,j,b,d);
                        }
                    
                }
                
            }
            
            mpi::all_reduce_function_tensor<mpi::op::sum>(kernel_1_ph, start);
            mpi::all_reduce_function_tensor<mpi::op::sum>(kernel_1_tph, start);
            mpi::all_reduce_function_tensor<mpi::op::sum>(kernel_1_pp, start);
            
            mpi::cout << "OK" << std::endl;
            
            //Construct Kernel-2 Functions
            mpi::cout << "Calculating Kernel-2 functions ... " << std::flush;
            
            std::vector<io::ctens> kernel_2_ph(nMatGB_kernel*nMatGF_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_2_pp(nMatGB_kernel*nMatGF_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            std::vector<io::ctens> kernel_2_tph(nMatGB_kernel*nMatGF_kernel, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            
            for (int om=start; om<end; om++){
                if (om >= nMatGB_kernel) break;
                for (int nu=0; nu<nMatGF_kernel; nu++){
                    
                    int const n = nu + om*nMatGF;
                    int const wf = greenOM.pos(omega_f_kernel(nu));
                    int const wb = greenOM.pos(omega_f_kernel(nu)-omega_b_kernel(om));
                    int const wc = greenOM.pos(omega_b_kernel(om)-omega_f_kernel(nu));
                    
                    int const n_susc_test = omega_b.pos(omega_b_kernel(om));
                    bool const do_conj = n_susc_test < 0 ? true : false;
                    int const om_susc = do_conj ? omega_b.pos(-omega_b_kernel(om)) : n_susc_test;
                    int const nu_neg = omega_f_kernel.pos(-omega_f_kernel(nu));
                    int const n_susc = do_conj ? nu_neg + om_susc*nMatGF : nu + om_susc*nMatGF;
                    
                    for(auto const& abcd : abcds){
                    
                        auto const a=abcd[1];
                        auto const b=abcd[2];
                        auto const c=abcd[3];
                        auto const d=abcd[4];
                        
                        std::string const entry = std::to_string(2*a)+"_"+std::to_string(2*b+1)+"_"+std::to_string(2*c)+"_"+std::to_string(2*d+1);
                        
                        kernel_2_ph[n].emplace(a,b,c,d,entry, -kernel_1_ph[om].at(a,b,c,d));
                        kernel_2_tph[n].emplace(a,b,c,d,entry, -kernel_1_tph[om].at(a,b,c,d));
                        kernel_2_pp[n].emplace(a,b,c,d,entry, -kernel_1_pp[om].at(a,b,c,d));
                        
                        for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                        for(std::size_t j = 0; j < jHybMatrix.size(); ++j){
                        
                            if (do_hph){
                                if (do_conj){
                                    kernel_2_ph[n](a,b,c,d)  -= 2.* std::conj(hedin_ph[n_susc].at(a,b,j,i))*U(i,c,j,d)/(green[wf](a,a)*green[wb](b,b));
                                } else {
                                    kernel_2_ph[n](a,b,c,d)  -= 2.* (hedin_ph[n_susc].at(a,b,j,i))*U(i,c,j,d)/(green[wf](a,a)*green[wb](b,b));
                                }
                            }
                                
                            
                            
                            if (do_htph){
                                if (do_conj) {
                                    kernel_2_tph[n](a,b,c,d) -= 2.*std::conj(hedin_tph[n_susc].at(a,i,j,d))*U(i,c,b,j)/(green[wf](a,a)*green[wb](d,d));
                                } else {
                                    kernel_2_tph[n](a,b,c,d) -= 2.*(hedin_tph[n_susc].at(a,i,j,d))*U(i,c,b,j)/(green[wf](a,a)*green[wb](d,d));
                                }
                            }
                                
                            
                            if (do_hpp){
                                if (do_conj){
                                    kernel_2_pp[n](a,b,c,d)  +=     std::conj(hedin_pp[n_susc].at(a,i,c,j))*U(j,i,b,d)/(green[wf](a,a)*green[wc](c,c));
                                } else {
                                    kernel_2_pp[n](a,b,c,d)  +=     (hedin_pp[n_susc].at(a,i,c,j))*U(j,i,b,d)/(green[wf](a,a)*green[wc](c,c));
                                }
                            }
                            
                            
                        }
                    }
                }
            }
            
            mpi::all_reduce_function_tensor<mpi::op::sum>(kernel_2_ph, start*nMatGF);
            mpi::all_reduce_function_tensor<mpi::op::sum>(kernel_2_tph, start*nMatGF);
            mpi::all_reduce_function_tensor<mpi::op::sum>(kernel_2_pp, start*nMatGF);
            
            mpi::cout << "OK" << std::endl;
            
            //Write results
            mpi::cout << "Outputting results ... " << std::flush;
            
            jsx::value jObservablesOut;
            jsx::value jKernel_0,jKernel_1,jKernel_2;
            
            jKernel_0["static"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_0);
            
            if (do_sph)
                jKernel_1["ph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_1_ph);
            if (do_stph)
                jKernel_1["tph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_1_tph);
            if (do_spp)
                jKernel_1["pp"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_1_pp);
            
            if (do_hph)
                jKernel_2["ph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_2_ph);
            if (do_htph)
                jKernel_2["tph"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_2_tph);
            if (do_hpp)
                jKernel_2["pp"] = func::write_functions<Value>(jParams, jHybMatrix, kernel_2_pp);
            
            jObservablesOut["kernel 0"] = std::move(jKernel_0);
            
            if (do_sph or do_spp or do_stph)
                jObservablesOut["kernel 1"] = std::move(jKernel_1);
            
            if (do_hph or do_hpp or do_htph)
                jObservablesOut["kernel 2"] = std::move(jKernel_2);
            
            mpi::cout << "OK" << std::endl;
            
            return jObservablesOut;
            
        }
        
        template<typename Value>
        jsx::value evaluateFullVertexFromKernels(jsx::value jParams, jsx::value const& jObservables) {
            
            auto const name = worm::Kernels::name;
            jsx::value jWorm = jParams(name);
            
            int b_cutoff =  jWorm.is("boson cutoff") ? std::max(1, int(jWorm("boson cutoff").int64())) : 1;
            int f_cutoff = jWorm.is("fermion cutoff") ? std::max(1, int(jWorm("fermion cutoff").int64())) : 50;
            int vb_cutoff =  jWorm.is("vertex boson cutoff") ? std::max(1, int(jWorm("vertex boson cutoff").int64())) : 1;
            int vf_cutoff = jWorm.is("vertex fermion cutoff") ? std::max(1, int(jWorm("vertex fermion cutoff").int64())) : 50;
            
            //Full vertex only for +/- range if complex
            int const pos_and_neg_freq = std::is_same<Value,double>::value ? 1 : 2;
            int const nMatGB = pos_and_neg_freq*vb_cutoff-pos_and_neg_freq+1;
            int const nMatGF = 2*vf_cutoff;
            int asymptotic_cutoff_l = jParams(name).is("asymptotic cutoff") ? jParams(name)("asymptotic cutoff").int64() : 10;
            int asymptotic_cutoff_l_4 = std::pow(asymptotic_cutoff_l,4);
            
            func::OmegaMap omega_f(nMatGF,false,true);
            func::OmegaMap omega_b(nMatGB,true,!std::is_same<Value,double>::value);
            
            //Kernels are always for full +/- frequency range
            int const nMatGB_kernel = 2*b_cutoff-1;
            int const nMatGF_kernel = 2*f_cutoff;
            
            func::OmegaMap omega_f_kernel(nMatGF_kernel,false,true);
            func::OmegaMap omega_b_kernel(nMatGB_kernel,true,true);
            
            jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
            std::vector<io::cmat> hyb; std::vector<io::Matrix<Value>> hybMoments;
            std::tie(hyb, hybMoments) = partition::func::get_hybridisation<Value>(jParams);
            
            //params::complete_impurity<Value>(jParams);
            imp::Tensor<Value> const U_tmp(jParams("hloc")("two body"),jHybMatrix.size());
            InteractionTensor<Value> const U(U_tmp,jHybMatrix.size());
            
            //Read in green functions and susc and hedin susceptibilities
            mpi::cout << "Reading in kernels ... " << std::flush;
            
            std::vector<io::ctens> kernel_1_ph = jObservables(worm::Kernels::name).is("kernel 1") and jObservables(worm::Kernels::name)("kernel 1").is("ph") ?
                meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 1")("ph"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size())
                : std::vector<io::ctens>(nMatGB_kernel);
            std::vector<io::ctens> kernel_1_tph = jObservables(worm::Kernels::name).is("kernel 1") and jObservables(worm::Kernels::name)("kernel 1").is("tph") ?
                meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 1")("tph"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size())
                : std::vector<io::ctens>(nMatGB_kernel);
            std::vector<io::ctens> kernel_1_pp = jObservables(worm::Kernels::name).is("kernel 1") and jObservables(worm::Kernels::name)("kernel 1").is("pp") ?
            meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 1")("pp"), jParams, jParams(cfg::susc_ph::Worm::name()), jHybMatrix, hyb.size())
                : std::vector<io::ctens>(nMatGB_kernel);;
            
            std::vector<io::ctens> kernel_2_ph = jObservables(worm::Kernels::name).is("kernel 2") and jObservables(worm::Kernels::name)("kernel 2").is("ph") ?
                meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 2")("ph"), jParams, jParams(worm::Kernels::name), jHybMatrix, hyb.size())
                : std::vector<io::ctens>(nMatGB_kernel*nMatGF_kernel);
            std::vector<io::ctens> kernel_2_tph = jObservables(worm::Kernels::name).is("kernel 2") and jObservables(worm::Kernels::name)("kernel 2").is("tph") ? meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 2")("tph"), jParams, jParams(worm::Kernels::name), jHybMatrix, hyb.size())
                : std::vector<io::ctens>(nMatGB_kernel*nMatGF_kernel);
            std::vector<io::ctens> kernel_2_pp = jObservables(worm::Kernels::name).is("kernel 2") and jObservables(worm::Kernels::name)("kernel 2").is("pp") ? meas::read_tensor_functions_from_obs<ut::complex>(jObservables(worm::Kernels::name)("kernel 2")("pp"), jParams, jParams(worm::Kernels::name), jHybMatrix, hyb.size())
                : std::vector<io::ctens>(nMatGB_kernel*nMatGF_kernel);
            
            mpi::cout << "OK" << std::endl;
            
            mpi::cout << "Reading in measured vertex ... " << std::flush;
            
            std::vector<io::ctens> measured_vertex = jParams.is(cfg::vertex_imprsum::Worm::name()) ?
            meas::read_tensor_functions_from_obs<ut::complex>(jObservables(cfg::vertex_imprsum::Worm::name())("full vertex"), jParams, jParams(cfg::vertex_imprsum::Worm::name()), jHybMatrix, hyb.size()) : std::vector<io::ctens>();
            
            auto const do_v = measured_vertex.size() > 0;
            if (!do_v){
                asymptotic_cutoff_l = 0;
                asymptotic_cutoff_l_4 = 0;
            }
            
            mpi::cout << "OK" << std::endl;
            
            //Construct Kernel-2 Functions
            mpi::cout << "Calculating Asymptotic Vertex from kernels ... " << std::flush;
            
            std::vector<io::ctens> asymptotic_vertex(nMatGB*nMatGF*nMatGF, io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));
            auto const abcds = do_v ? measured_vertex[0].ijkl() : construct_ijkls(jHybMatrix.size(), U);
            
            int const size = nMatGB;
            int const chunk = (size + mpi::number_of_workers() - 1)/mpi::number_of_workers();
            int const start = chunk*mpi::rank();
            int const end = chunk*(mpi::rank() + 1);
            
            for (int om=start; om<end; om++){
                if (om >= nMatGB) break;
                for (int nu1=0; nu1<nMatGF; nu1++)
                for (int nu2=0; nu2<nMatGF; nu2++){
                    
                    int const nu1_ph = omega_f_kernel.pos(omega_f(nu1));
                    int const nu1_tph = nu1_ph;
                    int const nu1_pp = nu1_ph;
                    
                    int const nu2_ph = omega_f_kernel.pos(omega_f(nu2));
                    int const nu2_tph = omega_f_kernel.pos(omega_f(nu1) - omega_b(om));
                    int const nu2_pp = nu2_ph;
                    
                    int const om_ph = omega_b_kernel.pos(omega_b(om));
                    int const om_tph = omega_b_kernel.pos(omega_f(nu1) - omega_f(nu2));
                    int const om_pp = omega_b_kernel.pos(omega_f(nu1) + omega_f(nu2) - omega_b(om));
                    
                    int const n = nu2 + nu1*nMatGF + om*nMatGF*nMatGF;
                    
                    int const n1_ph = nu1_ph + om_ph * nMatGF_kernel;
                    int const n1_tph = nu1_tph + om_tph * nMatGF_kernel;
                    int const n1_pp = nu1_pp + om_pp * nMatGF_kernel;
                    
                    int const n2_ph = nu2_ph + om_ph * nMatGF_kernel;
                    int const n2_tph = nu2_tph + om_tph * nMatGF_kernel;
                    int const n2_pp = nu2_pp + om_pp * nMatGF_kernel;
                                    
                    for(auto const& abcd :abcds){
                    
                        auto const a=abcd[1];
                        auto const b=abcd[2];
                        auto const c=abcd[3];
                        auto const d=abcd[4];
                        
                        if (do_v)
                            asymptotic_vertex[n].emplace(a,b,c,d,measured_vertex[0].entry(a,b,c,d), -2.*U(a,c,b,d));
                        else
                            asymptotic_vertex[n].emplace(a,b,c,d, std::to_string(2*a)+"_"+std::to_string(2*b+1)+"_"+std::to_string(2*c)+"_"+std::to_string(2*d+1) , -2.*U(a,c,b,d));
                            
                        
                        //PH
                        if (om_ph>=0){
                            asymptotic_vertex[n](a,b,c,d) += kernel_1_ph[om_ph].at(a,b,c,d);
                            if (nu1_tph>=0 and nu2_ph>=0){
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_ph[n1_ph].at(a,b,c,d) + kernel_2_ph[n2_ph].at(a,b,c,d);
                            }
                        }
                        
                        //TPH
                        if (om_tph>=0){
                            asymptotic_vertex[n](a,b,c,d) += kernel_1_tph[om_tph].at(a,b,c,d);
                            if (nu1_tph>=0 and nu2_tph>=0){
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_tph[n1_tph].at(a,b,c,d) + kernel_2_tph[n2_tph].at(a,b,c,d);
                            }
                            
                        }
                        
                        //PP
                        if (om_pp>=0){
                            asymptotic_vertex[n](a,b,c,d) += kernel_1_pp[om_pp].at(a,b,c,d);
                            if (nu1_pp>=0 and nu2_pp>=0){
                                asymptotic_vertex[n](a,b,c,d) += kernel_2_pp[n1_pp].at(a,b,c,d) + kernel_2_pp[n2_pp].at(a,b,c,d);
                            }
                            
                        }
                        
                    }
                
                }
            }
            
            mpi::all_reduce_function_tensor<mpi::op::sum>(asymptotic_vertex, start*nMatGF*nMatGF);
            
            mpi::cout << "OK" << std::endl;
            
            mpi::cout << "Combining asymptotic and measured vertices ... " << std::flush;
            
            auto const ref_worm =  do_v ? jParams(cfg::vertex_imprsum::Worm::name()) : jParams(worm::Kernels::name);
            auto const frequencies_meas =func::Frequencies<Value>(ref_worm);
            auto const& omega_f_meas = frequencies_meas.omega_f();
            auto const& omega_b_meas = frequencies_meas.omega_b();
            
            std::vector<io::ctens> combined_vertex(asymptotic_vertex.size(), io::ctens(jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size(), jHybMatrix.size()));

            for (int om=start; om<end; om++){
                if (om >= nMatGB) break;
                for (int nu1=0; nu1<nMatGF; nu1++)
                for (int nu2=0; nu2<nMatGF; nu2++){
                    
                    int const nu_prod = std::abs(omega_f(nu1) * (omega_f(nu1) - omega_b(om)) * (omega_f(nu2) - omega_b(om)) * omega_f(nu2));
                    int const is_static = 1;//!omega_b(om);
                    int const is_diag = nu1 == nu2;
                    int const use_asymptotic = nu_prod > asymptotic_cutoff_l_4*(is_static + is_diag - is_static*is_diag);
                    
                    int const om_meas = omega_b_meas.pos(omega_b(om));
                    int const nu1_meas = omega_f_meas.pos(omega_f(nu1));
                    int const nu2_meas = omega_f_meas.pos(omega_f(nu2));
                    
                    int const n = nu1 + nu2*nMatGF + om*nMatGF*nMatGF;
                    
                    int const n_meas = nu1_meas + nu2_meas*frequencies_meas.nMatGF() + om_meas*frequencies_meas.nMatGF()*frequencies_meas.nMatGF();
                    
                    for(auto const& abcd : abcds){
                    
                        auto const a=abcd[1];
                        auto const b=abcd[2];
                        auto const c=abcd[3];
                        auto const d=abcd[4];
                        
                        combined_vertex[n].emplace(a,b,c,d, asymptotic_vertex[0].entry(a,b,c,d), 0.);
                        
                        if (om_meas < 0 or nu1_meas < 0 or nu2_meas < 0 or use_asymptotic){
                            combined_vertex[n](a,b,c,d) = asymptotic_vertex[n](a,b,c,d);
                        } else {
                            combined_vertex[n](a,b,c,d) = measured_vertex[n_meas](a,b,c,d);
                        }
                            
                        
                    }
                    
                }
            }
            
            mpi::all_reduce_function_tensor<mpi::op::sum>(combined_vertex, start*nMatGF*nMatGF);
            
            mpi::cout << "OK" << std::endl;
            
            
            mpi::cout << "Outputting results ... " << std::flush;
            
            jsx::value jObservablesOut;
            
            jObservablesOut["full vertex"] = func::write_functions<Value>(jParams, jHybMatrix, combined_vertex);
            jObservablesOut["full vertex (asymptotic)"] = func::write_functions<Value>(jParams, jHybMatrix, asymptotic_vertex);
            
            mpi::cout << "OK" << std::endl;
            
            return jObservablesOut;
            
        }
        
        template jsx::value evaluateKernels<double>(jsx::value jParams, jsx::value const& jObservables);
        template jsx::value evaluateKernels<ut::complex>(jsx::value jParams, jsx::value const& jObservables);
        
        template jsx::value evaluateFullVertexFromKernels<double>(jsx::value jParams, jsx::value const& jObservables);
        template jsx::value evaluateFullVertexFromKernels<ut::complex>(jsx::value jParams, jsx::value const& jObservables);
    
        
        template <typename Value>
        std::vector<std::vector<int>> construct_ijkls(int const size, InteractionTensor<Value> const& U){
            std::vector<std::vector<int>> ijkls;
            int n=0;
            for (int i=0; i<size; i++)
                for (int j=0; j<size; j++)
                    for (int k=0; k<size; k++)
                        for (int l=0; l<size; l++)
                            if (std::abs(U(i,k,j,l)) > 1e-8 || (i==j and j==k and k==l) ){
                                ijkls.push_back(std::vector<int>({n,i,j,k,l}));
                                n+=1;
                            }
            
            return ijkls;
        };
        
        template std::vector<std::vector<int>> construct_ijkls<double>(int const size, InteractionTensor<double> const& U);
        template std::vector<std::vector<int>> construct_ijkls<ut::complex>(int const size, InteractionTensor<ut::complex> const& U);
        
    }
    
}





