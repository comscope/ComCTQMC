
#include "Algebra.h"
#include "MonteCarlo.h"

namespace mc {
    
    template<typename Mode, typename Value>
    void montecarlo(jsx::value jParams, jsx::value& jSimulation)
    {
        params::complete_impurity<Value>(jParams);
        
        data::Data<Value> data(jParams, Mode());
        data::setup_data<Mode>(jParams, data);
        
        obs::Observables<Value> observables;
        obs::setup_obs<Mode>(jParams, data, observables);
        
        mch::WangLandau<Value> wangLandau(jParams, data);
        
        std::vector<std::tuple<
        std::unique_ptr<imp::itf::Batcher<Value>>, //0
        std::unique_ptr<state::State<Value>     >, //1
        std::unique_ptr<mch::MarkovChain<Value> >, //2
        std::unique_ptr<mch::Scheduler          >  //3
        >> simulations;
        
        for(int stream = 0; stream < jSimulation.size(); ++stream) {
            
            simulations.emplace_back(
            std::unique_ptr<imp::itf::Batcher<Value>>(new imp::Batcher<Mode, Value>(8192)),
            std::unique_ptr<state::State<Value>     >(new state::State<Value>(jParams, data, jSimulation(stream)("config"), Mode())),
            std::unique_ptr<mch::MarkovChain<Value> >(new mch::MarkovChain<Value>(jParams, jSimulation(stream)("id").int64(), Mode())),
            std::unique_ptr<mch::Scheduler          >(new mch::Scheduler(false, mch::Phase::Initialize))
            );

            upd::setup_updates<Mode>(jParams, data, *std::get<1>(simulations.back()), *std::get<2>(simulations.back()), stream);
        }

        jSimulation["configs"] = jsx::array_t();
        
        meas::restart(jParams,jSimulation["measurements"]);
        
        std::int64_t thermSteps = 0, measSteps = 0, stream = 0;
        bool finalized = false; //if one stream gets has measured for x minutes / steps -- finalize all of the treams
        
        while(simulations.size()) {
            auto* batcher = std::get<0>(simulations[stream]).get();

            if(batcher->is_ready()) {
                auto& state       = std::get<1>(simulations[stream]);
                auto& markovChain = std::get<2>(simulations[stream]);
                auto& scheduler   = std::get<3>(simulations[stream]);

                try {
                    switch (scheduler->phase()) {
                        case mch::Phase::Step:
                            if(!markovChain->cycle(wangLandau, data, *state, *batcher)) break;
                                         
                            if(scheduler->done() or finalized) {
                                if(scheduler->thermalised() or finalized) {
                                    scheduler->phase() = mch::Phase::Finalize;
                                    finalized = true;
                                } else {
                                    if(jParams.is("measurement steps"))
                                        scheduler.reset(new mch::StepsScheduler(jParams("measurement steps").int64(), true, mch::Phase::Step));
                                    else
                                        scheduler.reset(new mch::TimeScheduler(jParams("measurement time").int64(), true, mch::Phase::Step, scheduler->overtime()));
                                    wangLandau.thermalised();
                                }
                                break;
                            }
                            
                            if(scheduler->thermalised()) {
                                ++measSteps;
                                
                                if(observables[state->worm().index()]->sample(data, *state))
                                    scheduler->phase() = mch::Phase::Sample;
                            } else
                                ++thermSteps;

                            break;
                            
                        case mch::Phase::Sample:
                            if(!observables[state->worm().index()]->cycle(data, *state, jSimulation["measurements"], *batcher)) break;
                            
                            scheduler->phase() = mch::Phase::Step;
                            
                            break;
                        
                        case mch::Phase::Finalize:
                            jSimulation["configs"].array().push_back(state->json());

                            simulations.erase(simulations.begin() + stream);
                            batcher = nullptr;

                            break;
                            
                        case mch::Phase::Initialize:
                            if(!markovChain->init(data, *state, *batcher)) break;

                            if(jParams.is("thermalisation steps"))
                                scheduler.reset(new mch::StepsScheduler(jParams("thermalisation steps").int64(), false, mch::Phase::Step));
                            else
                                scheduler.reset(new mch::TimeScheduler(jParams("thermalisation time").int64(), false, mch::Phase::Step));
                            
                            break;
                    }

                    if(batcher != nullptr) batcher->launch();
                }
                
                catch(ut::out_of_memory error) {
                    if(std::get<3>(simulations[stream])->phase() == mch::Phase::Sample)
                        throw std::runtime_error("MC: Fatal error, out of memory while sampling");
                    
                    std::cout << "MC: Markov Chain gets killed." << std::endl;

                    simulations.erase(simulations.begin() + stream);
                    
                    if(!simulations.size())
                        throw std::runtime_error("MC: all Markov Chain killed !");
                }
            }
            if(!(++stream < simulations.size())) stream = 0;
        }
        
        for(std::size_t space = 0; space < cfg::Worm::size(); ++space)
            if(observables[space] != nullptr)
                observables[space]->finalize(data, jSimulation["measurements"]);
        
        wangLandau.finalize(jSimulation["measurements"]);

        
        jSimulation["etas"] = wangLandau.etas();
        

        std::int64_t numberOfMarkovChains = jSimulation["configs"].size();
        
        mpi::reduce<mpi::op::sum>(numberOfMarkovChains, mpi::master);
        mpi::reduce<mpi::op::sum>(thermSteps,           mpi::master);
        mpi::reduce<mpi::op::sum>(measSteps,            mpi::master);

        jSimulation["info"] = jsx::object_t{
            { "number of mpi processes", mpi::number_of_workers() },
            { "number of markov chains", numberOfMarkovChains },
            { "thermalization steps",    thermSteps },
            { "measurement steps",       measSteps }
        };

    }
    
    template void montecarlo<imp::Host,double>(jsx::value jParams, jsx::value& jSimulation);
    template void montecarlo<imp::Host,ut::complex>(jsx::value jParams, jsx::value& jSimulation);
#ifdef MAKE_GPU_ENABLED
    template void montecarlo<imp::Device,double>(jsx::value jParams, jsx::value& jSimulation);
    template void montecarlo<imp::Device,ut::complex>(jsx::value jParams, jsx::value& jSimulation);
#endif
    
    
    template<typename Value>
    void statistics(jsx::value jParams, jsx::value& jSimulation) {
        if(mpi::number_of_workers() > 1 && jParams("error").string() != "none") {
            jsx::value jMeasurements = std::move(jSimulation("measurements"));
            
            meas::reduce(jSimulation("measurements"), jMeasurements, jSimulation("etas"), meas::All(), true);

            if(jParams("error").string() == "parallel") {
                
                //This very very slighly biases results towards 0 if it finds missing tensor elements -- ok for error only
                meas::check_missing_tensor_elements<Value>(jParams, jMeasurements);
                
                //measuring error accross workers means that each worker must calculate its own error
                //By default, we will not compute error of observables which are computed in parallel
                //i.e., those observables which might take a long time to compute serially
                jParams["serial evalsim"] = true;
                jParams["limited post-processing"] = !jParams("all errors").boolean();
                
                jsx::value jJackknifeError;
                
                meas::reduce(jJackknifeError, jMeasurements, jSimulation("etas"), meas::Jackknife(), false);
                jSimulation["error"] = evalsim::evalsim<Value>(jParams, jJackknifeError);
                meas::error(jSimulation("error"), meas::Jackknife());
                
                if (jParams.is("analytical continuation")){
                    

                    jParams["limited post-processing"] = false;
                    
                    jsx::value jBin;
                    meas::reduce(jBin, jMeasurements, jSimulation("etas"), meas::Rescale(), false);
                    jBin = evalsim::evalsim<Value>(jParams, jBin);
                    jBin = jBin["partition"]["aux green matsubara"];
                    auto jAvg = jBin;
                    auto jDif = jBin;
                    
                    meas::error(jAvg, meas::Average());
                    meas::subtract(jDif, jAvg);
                  
                    meas::error(jDif, meas::Covariance());
                    
                    jSimulation["covariance"] = jDif;
                    

                }
                
            } else if(jParams("error").string() == "serial") {
                
                meas::reduce(jMeasurements, jMeasurements, jSimulation("etas"), meas::Jackknife(), true);
                jSimulation["resample"] = std::move(jMeasurements);
                io::to_tagged_json(jSimulation("resample"));
                
            } else
                throw std::runtime_error("mc::statistics: invalid error option " + jParams("error").string());
            
        } else
            meas::reduce(jSimulation("measurements"), jSimulation("measurements"), jSimulation("etas"), meas::All(), true);
        
        io::to_tagged_json(jSimulation("measurements"));
    }
    
    template void statistics<double>(jsx::value jParams, jsx::value& jSimulation);
    template void statistics<ut::complex>(jsx::value jParams, jsx::value& jSimulation);
}




