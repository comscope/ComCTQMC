#ifndef CTQMC_INCLUDE_CONFIG_EXPANSION_H
#define CTQMC_INCLUDE_CONFIG_EXPANSION_H

#include <stdexcept>
#include <iostream>
#include <vector>

#include "../Utilities.h"
#include "../../../include/mpi/Utilities.h"
#include "../../../include/io/Vector.h"


namespace cfg {

    
    struct Entries : std::vector<ut::KeyType> {
        Entries() = default;
        Entries(Entries const&) = delete;
        Entries(Entries&&) = default; // ?????
        Entries& operator=(Entries const&) = delete;
        Entries& operator=(Entries&&) = default;
        ~Entries() = default;
        
        void insert(ut::KeyType key);
        void erase(ut::KeyType key);
    };
    
    
    struct Expansion : std::vector<Entries> {
        Expansion() = delete;
        Expansion(jsx::value jExpansion, std::size_t flavors);
        Expansion(Expansion const&) = delete;
        Expansion(Expansion&&) = delete;
        Expansion& operator=(Expansion const&) = delete;
        Expansion& operator=(Expansion&&) = delete;
        ~Expansion() = default;
        
        jsx::value json() const;
    };

}

#endif
