#ifndef CTQMC_INCLUDE_CONFIG_WORMS_H
#define CTQMC_INCLUDE_CONFIG_WORMS_H


#include <stdexcept>
#include <iostream>
#include <array>
#include <vector>

#include "Variant.h"
#include "Tuple.h"


namespace cfg {
    
    
    namespace partition {
        
        constexpr char name[] = "partition";
        
        using Worm = WormTuple<name>;
        
    }
    
    namespace green {
        
        constexpr char name[] = "green";
        
        using Worm = WormTuple<name, op, opDagg>;
        
    }
    
    namespace green_impr {
        
        constexpr char name[] = "green_impr";
        
        using Worm = WormTuple<name, opBulla, opDagg>;
        
    }
    
    namespace green_imprsum {
        
        constexpr char name[] = "green_imprsum";
        
        using Worm = WormTuple<name, opBullaSum, opDagg>;
        
    }
    
    namespace vertex {
        
        constexpr char name[] = "vertex";
        
        using Worm = WormTuple<name, op, opDagg, op, opDagg>;
        
    }
    
    namespace vertex_impr {
        
        constexpr char name[] = "vertex_impr";
        
        using Worm = WormTuple<name, opBulla, opDagg, op, opDagg>;
        
    }
    
    namespace vertex_imprsum {
        
        constexpr char name[] = "vertex_imprsum";
        
        using Worm = WormTuple<name, opBullaSum, opDagg, op, opDagg>;
        
    }
    
    namespace susc_ph {
        
        constexpr char name[] = "susc_ph";
        
        using Worm = WormTuple<name, bilinearPH, bilinearPH>;
        
    }
    
    namespace susc_pp {
        
        constexpr char name[] = "susc_pp";
        
        using Worm = WormTuple<name, bilinearHH, bilinearPP>;
        
    }
    
    namespace hedin_ph {
        
        constexpr char name[] = "hedin_ph";
        
        using Worm = WormTuple<name, op, opDagg, bilinearPH>;
        
    }
    
    namespace hedin_ph_impr {
        
        constexpr char name[] = "hedin_ph_impr";
        
        using Worm = WormTuple<name, opBulla, opDagg, bilinearPH>;
        
    }
    
    namespace hedin_ph_imprsum {
        
        constexpr char name[] = "hedin_ph_imprsum";
        
        using Worm = WormTuple<name, opBullaSum, opDagg, bilinearPH>;
        
    }
    
    namespace hedin_pp {
        
        constexpr char name[] = "hedin_pp";
        
        using Worm = WormTuple<name, op, op, bilinearPP>;
        
    }
    
    namespace hedin_pp_impr {
        
        constexpr char name[] = "hedin_pp_impr";
        
        using Worm = WormTuple<name, opBulla, op, bilinearPP>;
        
    }
    
    namespace hedin_pp_imprsum {
        
        constexpr char name[] = "hedin_pp_imprsum";
        
        using Worm = WormTuple<name, opBullaSum, op, bilinearPP>;
        
    }
    
    
    using Worm = WormVariant<
    
    partition::Worm,
    
    green::Worm,
    green_impr::Worm,
    green_imprsum::Worm,
    
    vertex::Worm,
    vertex_impr::Worm,
    vertex_imprsum::Worm,
    
    susc_ph::Worm,
    
    susc_pp::Worm,
    
    hedin_ph::Worm,
    hedin_ph_impr::Worm,
    hedin_ph_imprsum::Worm,
    
    hedin_pp::Worm,
    hedin_pp_impr::Worm,
    hedin_pp_imprsum::Worm
    
    >;
    
}

#endif
