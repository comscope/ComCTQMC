#ifndef CTQMC_INCLUDE_IMPURITY_PRODUCT_IMPL_H
#define CTQMC_INCLUDE_IMPURITY_PRODUCT_IMPL_H

#include "Product.h"

namespace imp {
    
    template<typename Mode, typename Value>
    template<typename... Args>
    typename Product<Mode,Value>::NodeType* Product<Mode,Value>::insert_impl(ut::KeyType const key, int const newHeight, Args&&... args) {
        if(!(0 <= key && key <= ut::KeyMax)) throw std::runtime_error("imp::Product::insert_impl: invalid key");

        std::vector<NodeType*> node(maxHeight_ + 1); NodeType* n = first_;
        std::vector<int> pos(maxHeight_ + 1); int p = 0;

        for(int l = std::min(height_, maxHeight_); l >= 0; --l) {
            for(; key > n->next[l]->key; n = n->next[l]) p += n->entries[l];
            node[l] = n; pos[l] = p; ++n->entries[l];
        }

        if(key == n->next[0]->key) {
            for(int l = std::min(height_, maxHeight_); l >= 0; --l) --node[l]->entries[l];
            return last_;
        }
        ++p; ++size_;

        for(int l = height_ + 1; l <= std::min(newHeight, maxHeight_); ++l) {
            node[l] = first_; pos[l] = 0; first_->next[l] = last_; first_->entries[l] = size_ + 1;
        }

        n = new NodeType(key, newHeight, std::forward<Args>(args)...);
        inserted_.push_back(n);

        for(int l = 0; l < newHeight; ++l) {
            n->next[l] = node[l]->next[l];
            node[l]->next[l] = n;

            int const diff = p - pos[l];
            n->entries[l] = node[l]->entries[l] - diff;
            node[l]->entries[l] = diff;
        }

        if(newHeight > height_) height_ = newHeight;
        
        if(node[0]->touch(0)) touched_.push_back(node[0]);
        for(int l = 1; l <= std::min(height_, maxHeight_); ++l)
            if(node[l - 1] != node[l])
                if(node[l]->touch(l)) touched_.push_back(node[l]);

        sign_ *= p%2 ? -1 : 1; return n;
    };

}

#endif


