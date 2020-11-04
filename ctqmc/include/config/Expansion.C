
#include "Expansion.h"

namespace cfg {

    
    void Entries::insert(ut::KeyType key) {
        std::vector<ut::KeyType>::insert(std::upper_bound(begin(), end(), key), key);
    };
    
    void Entries::erase(ut::KeyType key) {
        std::vector<ut::KeyType>::erase(std::lower_bound(begin(), end(), key));
    };
    
    Expansion::Expansion(jsx::value jExpansion, std::size_t flavors) {
        if(!jExpansion.is<jsx::null_t>()) {
            if(flavors != jExpansion.size())
                throw std::runtime_error("Expansion: invalid number of entries.");
            
            for(auto& jEntries : jExpansion.array()) {
                Entries entries;
                for(auto& key : jsx::at<io::ivec>(jEntries)) entries.push_back(key);
                push_back(std::move(entries));
            }
        } else
            resize(flavors);
    };
    
    
    jsx::value Expansion::json() const {
        jsx::array_t jExpansion;
        
        for(auto const& entries : *this) {
            io::ivec keys; keys.b64() = true;
            for(auto const& entry : entries) keys.push_back(entry);
            jExpansion.push_back(std::move(keys));
        }
        
        return jExpansion;
    };

}
