#include "Basic.h"

namespace cfg {
    
    bool operator<(Flavor const& lhs, Flavor const& rhs) {   // this guy is dangerous ....
        return lhs.flavor() < rhs.flavor();
    };
    
}
