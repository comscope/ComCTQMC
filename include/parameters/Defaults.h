#ifndef INCLUDE_PARAMETERS_STRUCTURES_H
#define INCLUDE_PARAMETERS_STRUCTURES_H

#include "../JsonX.h"

/*
Here we define default json dictionaries
*/
namespace params {

struct Defaults{
    
    jsx::value const& get() const {return defaults_;}
    
protected:
    jsx::value defaults_;
};

struct OneTimeWormDefaults : Defaults {

    OneTimeWormDefaults(){
        
        defaults_["cutoff"] = 50; // number of frequencies to measure
        defaults_["matsubara cutoff"] = 50; // number of frequencies to output in evalsim (used in legendre basis"
        defaults_["basis"] = "matsubara"; // basis in which to measure (only matsubara for susceptibilities atm)
        defaults_["meas"] = jsx::array_t({"imprsum"}); //use improved estimators/not (["imprsum"/""]) (improved estimators not implemented for susceptibilities)
        defaults_["sweep"] = 50;
        defaults_["store"] = 100;
        
    }
    
};

struct MultiTimeWormDefaults : Defaults {

    MultiTimeWormDefaults(){
        
        defaults_["fermion cutoff"] = 50; // number of frequencies to measure
        defaults_["boson cutoff"] = 10; // number of frequencies to measure
        defaults_["matsubara cutoff"] = 50; // number of frequencies to output in evalsim (used in legendre basis"
        defaults_["basis"] = "matsubara"; // basis in which to measure (only matsubara for susceptibilities atm)
        defaults_["meas"] = jsx::array_t({"imprsum"}); //use improved estimators/not (["imprsum"/""]) (improved estimators not implemented for susceptibilities)
        defaults_["sweep"] = 50;
        defaults_["store"] = 100;
        defaults_["full"] = true; //compute the full vertex -- only effects the "vertex x" worms and "kernels" block
        
    }
        
};


struct PartitionSpaceDefaults : Defaults {
    
    PartitionSpaceDefaults(){
        
        defaults_["green basis"] = "matsubara";
        defaults_["green bulla"] = true;
        defaults_["green matsubara cutoff"] = 10.001;
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

};

struct MainDefaults : Defaults{
    
    MainDefaults() = delete;
    MainDefaults(jsx::value const& jParams){
        defaults_["complex"] = jParams("hloc")("one body").is<jsx::object_t>();
        defaults_["restart"] = false;
        defaults_["partition fraction"] = 0.5;
        defaults_["sim per device"] = 0;
        defaults_["measurement time"] = 20;
        defaults_["thermalisation time"] = 5;
        defaults_["error"] = "parallel";
        defaults_["all error"] = false;
        defaults_["quad insert"] = false;
        defaults_["seed"] = 41085;
        defaults_["seed inc"] = 857;
        defaults_["expansion histogram"] = true;
        
    }
    
};

struct AllDefaults{
    
    AllDefaults(jsx::value const& jParams) : mainDefaults(jParams){
        
    }
    
public:
    MainDefaults mainDefaults;
    PartitionSpaceDefaults partitionSpaceDefaults;
    MultiTimeWormDefaults multiTimeWormDefaults;
    OneTimeWormDefaults oneTimeWormDefaults;
    
};

}

#endif
