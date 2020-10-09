#include  "../include/markovchain/MarkovChain.h"

namespace mch {
    
    std::int64_t select_seed(jsx::value const& jParams, std::int64_t const ID){
        
        auto inc = (jParams.is("seed increment") ? jParams("seed increment").int64() : 857)*ID;
        std::int64_t seed = 0;
        
        if (jParams.is("restart") and jParams("restart").boolean()){
            seed = static_cast<int64_t>(time(0));
            
        } else {
            seed = jParams.is("seed") ? jParams("seed").int64() : 41085;
            
        }
        
        return seed + inc;
    }
    
}


#include "../include/config/worm/Basic.h"

namespace cfg {

    bool operator<(Flavor const& lhs, Flavor const& rhs) {   // this guy is dangerous ....
        return lhs.flavor() < rhs.flavor();
    };

}

#include "../../include/options/Options.h"

namespace opt {
    
    Interaction get_interaction(jsx::value jBasis, jsx::value jTwoBody) {
        if(jBasis("orbitals").is<jsx::string_t>()) {
            
            sphericalharmonics::Basis basis(jBasis);
            
            if(jTwoBody("parametrisation").string() == "slater-condon")
                return Interaction(sphericalharmonics::SlaterCondon(basis, jTwoBody));
            else
                throw std::runtime_error("opt: parametrisation" + jTwoBody("parametrisation").string() + " not defined for spherical harmonics basis");
            
        } else if(jBasis("orbitals").is<jsx::int64_t>()) {
            
            model::Basis basis(jBasis);
            
            if(jTwoBody("parametrisation").string() == "kanamori")
                return Interaction(model::Kanamori(basis, jTwoBody));
            else
                throw std::runtime_error("opt: parametrisation" + jTwoBody("parametrisation").string() + " not defined for model basis");
            
        } else
            throw std::runtime_error("opt: orbitals not defined");
    }
    
    Observable get_observable(jsx::value const& jBasis, std::string const name) {
        if(jBasis("orbitals").is<jsx::string_t>()) {

            if(jBasis("type").string() == "product real" || jBasis("type").string() == "product") {
                if(name == "S2")
                    return Observable(sphericalharmonics::S2(sphericalharmonics::Basis(jBasis)));
            } else if(jBasis("type").string() == "product imag") {
                if(name == "S2")
                    return Observable(sphericalharmonics::S2(sphericalharmonics::Basis(jBasis)));
                //else if(obs.first == "L2")
                //    return Observable(sphericalharmonics::L2(l), transformation);
            } else if(jBasis("type").string() == "coupled") {
                if(name == "J2")
                    return Observable(sphericalharmonics::J2(sphericalharmonics::Basis(jBasis)));
            }
            
            throw std::runtime_error("opt: observable " + name + " not defined for spherical harmonics");

        } else if(jBasis("orbitals").is<jsx::int64_t>()) {

            if(name == "S2")
                return Observable(model::S2(model::Basis(jBasis)));
            
            throw std::runtime_error("opt: observable " + name + " not defined for model basis");
            
        } else
            
            throw std::runtime_error("opt: orbitals not defined");
    }
    
}

#include "../../include/measurements/Measurements.h"
namespace meas {
    
    std::int64_t reduce_steps(std::int64_t steps, All) {
        mpi::reduce<mpi::op::sum>(steps, mpi::master); return steps;
    };
    
    std::int64_t reduce_steps(std::int64_t steps, Jackknife) {
        auto temp = steps; mpi::all_reduce<mpi::op::sum>(temp); return temp - steps;
    };
    
    std::vector<std::string> split_by_char(std::string const& string, char const c){
        std::stringstream temp(string);
        std::string segment;
        std::vector<std::string> seglist;

        while(std::getline(temp, segment, c))
        {
           seglist.push_back(segment);
        }
        
        return seglist;
    }
        
    
    void restart(jsx::value const& jIn, jsx::value & jMeas, double const eta, int const samples){
        
        if(jIn.is<jsx::array_t>()) {
            if(!(jMeas.is<jsx::array_t>() && jMeas.size() == jIn.size())) jMeas = jsx::array_t(jIn.size());
            int index = 0; for(auto const& jEntry : jIn.array()) restart(jEntry, jMeas[index++], eta, samples);
            
        } else if (jIn.is("io::cvec")){ //for some reason, jsx::value.is<io::rvec>() doesn't recognize the io::cvec?)
            auto jData = jIn("io::cvec");
            std::vector<ut::complex> data = jsx::at<io::cvec>(jData);
            for (auto& x : data) x*=samples*eta;
            jMeas << meas::fix(data, samples);
            
        } else if (jIn.is("io::rvec")) { //jIn.is<io::rvec>()
            auto jData = jIn("io::rvec");
            std::vector<double> data = jsx::at<io::rvec>(jData);
            for (auto& x : data) x*=samples*eta;
            jMeas << meas::fix(data, samples);
            
        } else if (jIn.is<jsx::object_t>()){
            for (auto const& jEntry : jIn.object()){
                if (jEntry.first != "steps" and jEntry.first != "eta"){ //not real measurements -- WangLandau deals with them
                    if (jEntry.first != "expansion histogram") //the only vecvar I've encountered -- TODO: an automatic way to do this
                        restart(jEntry.second, jMeas[jEntry.first], eta, samples);
                }
            }
        }
        
    }

    void restart(jsx::value & jParams, jsx::value & jMeasurements){
        
        if(jParams.is("restart") and jParams("restart").boolean()){
            
            mpi::cout << "Initializing measurements from previous run" << std::endl;
            auto const names = cfg::Worm::get_names();
            
            //The measurements stored during ctqmc are not yet normalized by the relative size of the configuration spaces
            //or the average sign, but the reported measurements in params.meas.json are
            //so, we must go back to the raw measurement
            //this is the 1st part of refined -> raw (see reduce, e.g.)
            auto const pName = cfg::partition::Worm::name();
            auto jSign = jParams("measurements")(pName)("sign")("io::rvec");
            std::vector<double> sign = jsx::at<io::rvec>(jSign);
            auto const pEta = jParams("measurements")(pName)("eta").real64();
            auto const pSteps = jParams("measurements")(pName)("steps").int64();
            auto const signxZp = sign.at(0)*pSteps/pEta;
            
            for (auto & jIn : jParams("measurements").object()){ //loop through configuration spaces
                
                //2nd part of the refined -> raw
                auto const Zw = jIn.second("steps").int64()/jIn.second("eta").real64();
                
                //weight previous measurements of by # of samples taken
                std::int64_t sweep = 0;
                if (jIn.first == cfg::partition::Worm::name()){
                    std::int64_t const sweepB = jParams(jIn.first).is("sweepB") ? jParams(jIn.first)("sweepB").int64() : 250;
                    std::int64_t const sweepA = jParams(jIn.first).is("sweepA") ? jParams(jIn.first)("sweepA").int64() : 50;
                    sweep = std::max(sweepB, sweepA);
                    //This will slow convergence for the smaller of sweepB/sweepB (by default, the A-type partition obs)
                    //but this is better than preventing convergence for the other type.
                    //Consider using sweepB only if susceptibility bulla or green bulla set
                    //TODO: It would be better if params.meas.json output the sweep / var / fix of each observable so that restart could use it
                } else {
                    sweep = jParams(jIn.first).is("sweep") ? jParams(jIn.first)("sweep").int64() : 50;
                }
                
                std::size_t samples = jIn.second("steps").int64() / (mpi::number_of_workers() * sweep);
                
                restart(jIn.second, jMeasurements[jIn.first], signxZp/Zw, samples);
                
            }
            
            jParams["measurements"]=jsx::empty_t(); //free memory -- old measurements
            
            mpi::cout << "Done initializing measurements from previous run" << std::endl;
        }
    }
    
}

#include "../../include/atomic/Generate.h"
namespace ga {
    
    FlavorState::FlavorState() : sign_(1) {};
    
    FlavorState::FlavorState(FlavorState const& state) : std::bitset<8*sizeof(State)>(state), sign_(state.sign()) {};
    
    State FlavorState::state() const { return this->to_ullong();};
    int FlavorState::sign() const { return sign_;};
    int& FlavorState::sign() { return sign_;};
    
    FlavorState psi(int flavor, FlavorState const& state) {
        FlavorState temp = state;
        if(state.test(flavor)) {
            temp[flavor] = 0; if((state >> (flavor + 1)).count()%2) temp.sign() *= -1;
        } else
            temp.sign() = 0;
        return temp;
    };
    
    FlavorState psiDagg(int flavor, FlavorState const& state) {
        FlavorState temp = state;
        if(state.test(flavor))
            temp.sign() = 0;
        else {
            temp[flavor] = 1; if((state >> (flavor + 1)).count()%2) temp.sign() *= -1;
        }
        return temp;
    };
    
    
    Join::Join(std::size_t size) : labels_(size) {
        std::iota(labels_.begin(), labels_.end(), 0);
    };
    void Join::join(State state1, State state2) {
        state1 = find_representative(state1);
        state2 = find_representative(state2);
        state1 < state2 ? labels_[state2] = state1 : labels_[state1] = state2;
    };
    std::size_t Join::clean() {
        for(auto & label : labels_) label = find_representative(label);
        std::map<State, State> new_label;
        for(auto & label : labels_) {
            if(!new_label.count(label)) {
               State temp = new_label.size();
               new_label[label] = temp;
            }
            label = new_label[label];
        }
        return new_label.size();
    };
    State Join::label(State state) {
        return labels_[state];
    };
    State Join::find_representative(State label) {
        while(label != labels_[label]) label = labels_[label];
        return label;
    };
    
    io::rvec get_sector_qn(BlockStates const& blockStates,
                             std::vector<double> const& qn,
                           bool const throw_error)
    {
        io::rvec Qn;
        double const eps = 1E-4;
        
        for(auto const& states : blockStates) {
            double value = .0;
            
            for(std::size_t f = 0; f < qn.size(); ++f)
                if(FlavorState(states.front()).test(f)) value += qn[f];
            
            for(auto const state : states) {
                double temp = .0;
                
                for(std::size_t f = 0; f < qn.size(); ++f)
                    if(FlavorState(state).test(f)) temp += qn[f];
                
                if(std::abs(value  - temp) > eps){
                    if (throw_error) throw std::runtime_error("Problem with quantum numbers.");
                    else {mpi::cout << " not a good qn; removing from list " ; return io::rvec();}
                }
                    
            }
            
            Qn.push_back(value);
        }
        
        return Qn;
    };
    
    
    
    BlockStates get_block_states(jsx::value const& jHloc)
    {
        BlockStates blockStates;
        
        for(auto const& jStates : jHloc("block states").array()) {
            States states;
            for(auto const& jState : jStates.array())
                states.push_back(jState.int64());
            blockStates.push_back(std::move(states));
        }
        
        return blockStates;
    };
    
    io::rvec construct_sector_qn(jsx::value const& jHloc, jsx::value jqn, bool throw_error)
    {
        return get_sector_qn(get_block_states(jHloc), jsx::at<io::rvec>(jqn), throw_error);
    };
    
}

#include "../include/Params.h"
namespace params{
    
    void complete_worm(jsx::value& jParams, std::string const worm)
    {
        if(jParams.is(worm)) {
            jsx::value jEntry = jParams(worm);
            jParams.object().erase(worm);
            
            jsx::value jMeas = jEntry("meas");
            jEntry.object().erase("meas");
            
            std::map<std::string, std::string> map{{"",""}, {"impr", " impr"}, {"imprsum", " imprsum"}};
            
            for(auto& meas : jMeas.array()) {
                if(!map.count(meas.string()))
                    throw std::runtime_error("params::complete_worms: invalid meas option " + meas.string());
                meas = map.at(meas.string());
            }
            
            for(auto meas : jMeas.array()) {
                jParams[worm + meas.string()] = jEntry;
                if(meas.string() != "" && jParams.is("dyn"))
                    jParams[worm] = jEntry;
            }
            
            for(auto meas : jMeas.array()) {
                jParams(worm + meas.string())["static"] = jsx::null_t();
                if(meas.string() != "" && jParams.is("dyn"))
                    jParams(worm)["dynamic"] = jsx::null_t();
            }
        }
    }
    
    void complete_worms(jsx::value& jParams)
    {
        complete_worm(jParams, cfg::green::name);
        
        complete_worm(jParams, cfg::vertex::name);
        
        complete_worm(jParams, cfg::hedin_ph::name);
        
        complete_worm(jParams, cfg::hedin_pp::name);
        
        if(jParams.is(cfg::susc_ph::name))
           jParams(cfg::susc_ph::name)["static"] = jsx::null_t();
        
        if(jParams.is(cfg::susc_pp::name))
            jParams(cfg::susc_pp::name)["static"] = jsx::null_t();
    };
    
}
