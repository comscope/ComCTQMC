
#include "Node.h"
#include "../Algebra.h"

namespace imp {
    
    //Scheisse das muess verbessert werde chunt ja chei sau druus da ...
    
    template<typename NodeValue>
    Node<NodeValue>::~Node() {
        delete[] entries; delete[] next;
    };
    
    template<typename NodeValue>
    void Node<NodeValue>::accept() {
        std::memcpy(next + height + touched, next + touched, (height - touched)*sizeof(Node*));
        std::memcpy(entries + height + touched, entries + touched, (height - touched)*sizeof(int));
        value.accept(); touched = height;
    };
    
    template<typename NodeValue>
    void Node<NodeValue>::reject() {
        std::memcpy(next + touched, next + height + touched, (height - touched)*sizeof(Node*));
        std::memcpy(entries + touched, entries + height + touched, (height - touched)*sizeof(int));
        value.reject(); touched = height;
    };
    
    template<typename NodeValue>
    bool Node<NodeValue>::touch(int level) {
        bool temp = (touched == height);
        touched = std::min(level, touched);
        return temp;
    };
    
    
    template<typename Mode, typename Value>
    NodeValue<Mode,Value>::NodeValue(Access<NodeValue const> node, itf::Operator<Value> const* op0, int flavor, EigenValues<Mode> const& eig) :
    op0(&get<Mode, Value>(*op0)), flavor(flavor),
    eig_(eig), node_(node),
    prop_(nullptr), propTry_(nullptr),
    ops_(new Operator<Mode, Value>*[2*node_.height()]), opsTry_(ops_ + node_.height()) {
        for(int l = 0; l < 2*node_.height(); ++l) ops_[l] = nullptr;
    };
    
    template<typename Mode, typename Value>
    NodeValue<Mode,Value>::~NodeValue() {
        for(int l = 0; l < 2*node_.height(); ++l) delete ops_[l];
        delete[] ops_; delete propTry_; delete prop_;
    };
    
    template<typename Mode, typename Value>
    Propagator<Mode> const* NodeValue<Mode,Value>::prop() const {
        auto prop = node_.touched() ? prop_ : propTry_;
        if(prop == nullptr) throw std::runtime_error("imp::Value::prop: null pointer");
        return prop;
    };
    
    template<typename Mode, typename Value>
    Operator<Mode, Value> const* NodeValue<Mode,Value>::op(int l) const {
        auto op = l < node_.touched() ? ops_[l] : opsTry_[l];
        if(op == nullptr) throw std::runtime_error("imp::Value::op: null pointer");
        return op;
    };
    
    template<typename Mode, typename Value>
    std::unique_ptr<Operator<Mode, Value>> NodeValue<Mode,Value>::get_op(int l) const {
        auto& op = l < node_.touched() ? ops_[l] : opsTry_[l];
        if(op == nullptr) throw std::runtime_error("imp::Value::get_op: null pointer");
        auto temp = op; op = nullptr; return std::unique_ptr<Operator<Mode, Value>>(temp);
    };
    
    template<typename Mode, typename Value>
    Propagator<Mode>* NodeValue<Mode,Value>::prop() {
        auto& prop = node_.touched() ? prop_ : propTry_;
        return prop ? prop : prop = new Propagator<Mode>(-(node_.next(0).key() - node_.key())*ut::beta()/ut::KeyMax, eig_);  //promotion stuff ...
    };
    
    template<typename Mode, typename Value>
    Operator<Mode, Value>* NodeValue<Mode,Value>::op(int l) {
        auto& op = l < node_.touched() ? ops_[l] : opsTry_[l];
        return op ? op : op = new Operator<Mode, Value>(eig_);
    };
    
    template<typename Mode, typename Value>
    void NodeValue<Mode,Value>::accept() {
        if(!node_.touched()) {
            delete prop_; prop_ = propTry_; propTry_ = nullptr;
        }
        for(int l = std::max(node_.touched(), 1); l < node_.height(); ++l) {
            delete ops_[l]; ops_[l] = opsTry_[l]; opsTry_[l] = nullptr;
        }
    };
    
    template<typename Mode, typename Value>
    void NodeValue<Mode,Value>::reject() {
        if(!node_.touched()) {
            delete propTry_; propTry_ = nullptr;
        }
        for(int l = std::max(node_.touched(), 1); l < node_.height(); ++l) {
            delete opsTry_[l]; opsTry_[l] = nullptr;
        }
    };
    
    template struct NodeValue<imp::Host,double>;
    template struct NodeValue<imp::Host,ut::complex>;
    
    template struct Node<NodeValue<imp::Host,double>>;
    template struct Node<NodeValue<imp::Host,ut::complex>>;
    
#ifdef MAKE_GPU_ENABLED
    template struct NodeValue<imp::Device,double>;
    template struct NodeValue<imp::Device,ut::complex>;
    
    template struct Node<NodeValue<imp::Device,double>>;
    template struct Node<NodeValue<imp::Device,ut::complex>>;
#endif
    
}

