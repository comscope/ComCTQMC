#ifndef CTQMC_INCLUDE_IMPURITY_NODE_IMPL_H
#define CTQMC_INCLUDE_IMPURITY_NODE_IMPL_H

#include "Node.h"

namespace imp {

    template<typename NodeValue>
    template<typename... Args>
    Node<NodeValue>::Node(ut::KeyType key, int height, Args&&... args) :
        key(key), height(height),
        next(new Node*[2*height]),
        entries(new int[2*height]),
        touched(0),
        value(this, std::forward<Args>(args)...) { };
        

}

#endif


