#ifndef CTQMC_INCLUDE_IMPURITY_NODE_H
#define CTQMC_INCLUDE_IMPURITY_NODE_H

#include <iostream>
#include <random>

#include "Diagonal.h"
#include "Operators.h"
#include "../Utilities.h"

namespace imp {
    
    //Scheisse das muess verbessert werde chunt ja chei sau druus da ...

    template<typename NodeValue>
    struct Node {
        template<typename... Args>
        Node(ut::KeyType key, int height, Args&&... args);
        ~Node();
        
        ut::KeyType const key; int const height;
        Node** const next; int* const entries;
        
        void accept();
        void reject();
        bool touch(int level);

        int touched; NodeValue value;
    };
    
    
    template<typename NodeValue, typename F = void>
    struct Access {
        typedef Node<typename std::remove_cv<NodeValue>::type> NodeType;
        
        Access() = default;
        Access(NodeType* node) : node_(node) {};
        Access(Access const&) = default;
        Access(Access&&) = default;
        Access& operator=(Access const&) = default;
        Access& operator=(Access&&) = default;
        ~Access() = default;
        
        int height() const { return node_->height;};
        ut::KeyType key() const { return node_->key;};
        Access next(int level) const { return node_->next[level];};
        int entries(int level) const { return node_->entries[level];};
        int touched() const { return node_->touched;};
        
        NodeValue* operator->() const { return &node_->value;};
        
    private:
        NodeType* node_;
        
        friend inline bool operator==(Access<NodeValue, F> const& lhs, Access<NodeValue, F> const& rhs) { return lhs.node_ == rhs.node_;};
        friend inline bool operator!=(Access<NodeValue, F> const& lhs, Access<NodeValue, F> const& rhs) { return lhs.node_ != rhs.node_;};
        
        friend F;
    };
    
    
    template<typename Mode, typename Value>
    struct NodeValue {
        NodeValue(Access<NodeValue const> node, itf::Operator<Value> const* op0, int flavor, EigenValues<Mode> const& eig);
        ~NodeValue();
        
        Operator<Mode, Value> const* const op0;
        int const flavor; 
        
        
        Propagator<Mode> const* prop() const;
        Operator<Mode, Value> const* op(int l) const;
        std::unique_ptr<Operator<Mode, Value>> get_op(int l) const;

        Propagator<Mode>* prop();
        Operator<Mode, Value>* op(int l);
        
        void accept();
        void reject();
        
    private:
        EigenValues<Mode> const& eig_;
        Access<NodeValue const> node_;
        
        Propagator<Mode>* prop_; Propagator<Mode>* propTry_;
        Operator<Mode, Value>** const ops_; Operator<Mode, Value>** const opsTry_;
    };
}

#include "Node.impl.h"

#endif


