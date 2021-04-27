
#include "Defaults.h"

namespace params {

    jsx::value const& Defaults::get() const {return defaults_;}

    OneTimeWormDefaults::OneTimeWormDefaults(){
        
        defaults_["cutoff"] = 50; // number of frequencies to measure
        defaults_["matsubara cutoff"] = 50; // number of frequencies to output in evalsim (used in legendre basis"
        defaults_["basis"] = "matsubara"; // basis in which to measure (only matsubara for susceptibilities atm)
        defaults_["meas"] = jsx::array_t({"imprsum"}); //use improved estimators/not (["imprsum"/""]) (improved estimators not implemented for susceptibilities)
        defaults_["sweep"] = 50;
        defaults_["store"] = 100;
        defaults_["diagonal"] = false;
        
    }
    
    MultiTimeWormDefaults::MultiTimeWormDefaults(){
        
        defaults_["fermion cutoff"] = 50; // number of frequencies to measure
        defaults_["boson cutoff"] = 10; // number of frequencies to measure
        defaults_["matsubara cutoff"] = 50; // number of frequencies to output in evalsim (used in legendre basis"
        defaults_["basis"] = "matsubara"; // basis in which to measure (only matsubara for susceptibilities atm)
        defaults_["meas"] = jsx::array_t({"imprsum"}); //use improved estimators/not (["imprsum"/""]) (improved estimators not implemented for susceptibilities)
        defaults_["sweep"] = 50;
        defaults_["store"] = 100;
        defaults_["full"] = true; //compute the full vertex -- only effects the "vertex x" worms and "kernels" block
        defaults_["diagonal"] = false;
        
    }

    PartitionSpaceDefaults::PartitionSpaceDefaults(){
        
        defaults_["green basis"] = "matsubara";
        defaults_["green bulla"] = true;
        defaults_["green matsubara cutoff"] = 10;
        defaults_["green legendre cutoff"] = 100;
        
        defaults_["density matrix precise"] = false;
        
        defaults_["occupation susceptibility"] = false;
        defaults_["occupation susceptibility direct"] = false;
        defaults_["occupation susceptibility bulla"] = false;
        defaults_["quantum number susceptibility"] = false;
        defaults_["susceptibility cutoff"] = 10;
        defaults_["susceptibility tail"] = 100;
        
        defaults_["probabilities"] = jsx::array_t({"energy","N"});
        
        defaults_["sweepA"] = 50;
        defaults_["storeA"] = 100;
        
        defaults_["sweepB"] = 250;
        defaults_["storeB"] = 20;
        
    }

    MainDefaults::MainDefaults(jsx::value const& jParams){
        
        defaults_["restart"] = false;
        defaults_["partition fraction"] = 0.5;
        defaults_["sim per device"] = 0;
        defaults_["measurement time"] = 20;
        defaults_["thermalisation time"] = 5;
        defaults_["error"] = "parallel";
        defaults_["all errors"] = false;
        defaults_["quad insert"] = false;
        defaults_["seed"] = 41085;
        defaults_["seed inc"] = 857;
        defaults_["expansion histogram"] = true;
        
    }
    

    MainRequired::MainRequired(){
        defaults_["complex"] = jsx::empty_t();
        defaults_["beta"] = jsx::empty_t();
        defaults_["mu"] = jsx::empty_t();
        defaults_[cfg::partition::Worm::name()] = jsx::empty_t();
        defaults_["hloc"] = jsx::array_t({"one body", "two body"});
        defaults_["hybridisation"] = jsx::array_t({"functions", "matrix"});
    }
    
    void MainRequired::init_types(){
        defaults_["complex"] = false;
        defaults_["beta"] = 0.;
        defaults_["mu"] = 0.;
        
        defaults_[cfg::partition::Worm::name()] = jsx::object_t();
        
        //leave the two body check for the hloc generation, which is both more complicated and also throws reasonable errors
        if (defaults_["complex"].boolean()){
            defaults_["hloc"]["one body"] = io::PrettyMatrix<ut::complex>();
        } else {
            defaults_["hloc"]["one body"] = io::PrettyMatrix<double>();
        }
        
            
        defaults_["hybridisation"]["functions"] = "hyb.json";
        defaults_["hybridisation"]["matrix"] = jsx::array_t({ jsx::array_t({""}) });
    }


    AllDefaults::AllDefaults(jsx::value const& jParams) : mainDefaults(jParams){
        
         listOfDefaultWorms = std::vector<std::string>({
             cfg::susc_ph::Worm::name(),
             cfg::susc_pp::Worm::name(),
             cfg::hedin_ph::Worm::name(),
             cfg::hedin_pp::Worm::name(),
             cfg::vertex::Worm::name()
         });
    }
    
}

