#ifndef INCLUDE_PARAMETERS_DEFAULTS_H
#define INCLUDE_PARAMETERS_DEFAULTS_H

#include "../JsonX.h"
#include "../../ctqmc/include/config/Worms.h"

/*
Here we define default json dictionaries for most input blocks
Some blocks, like the basis and interaction, are handled in the options namespace
*/
namespace params {

    struct Defaults{
        
        jsx::value const& get() const;
        
    protected:
        jsx::value defaults_;
    };

    struct OneTimeWormDefaults : Defaults {
        
        OneTimeWormDefaults();
        
    };

    struct MultiTimeWormDefaults : Defaults {
        
        MultiTimeWormDefaults();
        
    };


    struct PartitionSpaceDefaults : Defaults {
        
        PartitionSpaceDefaults();
        
    };

    struct MainDefaults : Defaults{
        
        MainDefaults() = delete;
        MainDefaults(jsx::value const& jParams);
        
    };


    struct MainRequired : Defaults{
        
        //To make sure these things exist at all
        MainRequired();
        
        //To make sure they are the right datatypes
        void init_types();
        
    };

    struct AllDefaults{
        
        AllDefaults(jsx::value const& jParams);
        
    public:
        
        std::vector<std::string> listOfDefaultWorms;
        
        MainDefaults mainDefaults;
        PartitionSpaceDefaults partitionSpaceDefaults;
        MultiTimeWormDefaults multiTimeWormDefaults;
        OneTimeWormDefaults oneTimeWormDefaults;
        
    };

}

#include "Defaults.impl.h"

#endif
