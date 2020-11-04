#ifndef CTQMC_INCLUDE_IMPURITY_PRODUCT_H
#define CTQMC_INCLUDE_IMPURITY_PRODUCT_H

#include <iostream>
#include <random>

#include "Algebra.h"
#include "Node.h"
#include "../Utilities.h"
#include "../../../include/JsonX.h"

namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct Product {
            virtual int size() const = 0;
            virtual int height() const = 0;
            virtual bool insert(ut::KeyType const key, int flavor) = 0;
            virtual void erase(ut::KeyType key) = 0;
            virtual int accept() = 0;
            virtual void reject() = 0;
            virtual ~Product() = default;
        };
        
    };
    
    
    // Todo: remove flavor from node 
    
    template<typename Mode, typename Value>
    struct Product : itf::Product<Value> {
        using NodeType  = Node<NodeValue<Mode, Value>>;
        
        using AccessType  = Access<NodeValue<Mode, Value>, Product>;
        using CAccessType = Access<NodeValue<Mode, Value> const, Product>;
        
        Product() = delete;
        Product(jsx::value const& jParams, itf::EigenValues const& eig, itf::Operator<Value> const& ide, itf::Operators<Value> const& ops);
        Product(Product const&) = delete;
        Product(Product&&) = delete;
        Product& operator=(Product const&) = delete;
        Product& operator=(Product&&) = delete;
        ~Product();

        int size() const { return size_;};
        int height() const { return height_;};

        CAccessType first() const { return first_;};
        CAccessType last() const { return last_;};
        
        bool insert(ut::KeyType const key, int flavor);
        CAccessType insert(ut::KeyType const key, itf::Operator<Value> const* op, int flavor);
        CAccessType insert(ut::KeyType const key, itf::Operator<Value> const* op, int flavor, int height);
        
        void erase(ut::KeyType key);
        
        std::vector<int> map(CAccessType cbegin, int level, std::vector<int> const& sectors);
        void multiply(CAccessType cbegin, int level, int sec, itf::Batcher<Value>& batcher);
		
        int accept();
        void reject();
        
        Matrix<Mode, Value>& bufferA() { return *bufferA_;};
        Matrix<Mode, Value>& bufferB() { return *bufferB_;};
        Matrix<Mode, Value>& bufferC() { return *bufferC_;};

    private:
        EigenValues<Mode> const& eig_;
        Operator<Mode, Value> const& ide_;
        Operators<Mode, Value> const& ops_;
        
		ut::RandomNumberGenerator<std::mt19937, std::uniform_real_distribution<double> > urng_;
		
		double const prob_; 
		double const baseProb_;
		int const maxHeight_;	
		int height_, heightBackup_;
		int size_, sizeBackup_;
        int sign_;

		NodeType* const first_; NodeType* const last_;		

		std::vector<NodeType*> touched_, inserted_, erased_;
        
        SectorNormPtrs all_;
        
        std::vector<Matrix<Mode, Value> const*> ptrMat_;
        std::vector<Vector<Mode> const*> ptrProp_;
        std::unique_ptr<Matrix<Mode, Value>> bufferA_, bufferB_, bufferC_;

        int random_height();
        
        template<typename... Args>
        NodeType* insert_impl(ut::KeyType const key, int const newHeight, Args&&... args);
        
        void erase_impl(ut::KeyType const key);
        
        static void map_impl(AccessType it, int level, SectorNormPtrs& requested);
        
        static void multiply_impl(AccessType const begin, int const level, int const sec, Matrix<Mode, Value> const** const mat, Vector<Mode> const** const prop, Matrix<Mode, Value>* A, Matrix<Mode, Value>* B, itf::Batcher<Value>& batcher);
    };
    
    template<typename Mode, typename Value> Product<Mode, Value>& get(itf::Product<Value>& productItf) {
        return static_cast<Product<Mode, Value>&>(productItf);
    };
    
    template<typename Mode, typename Value> Product<Mode, Value> const& get(itf::Product<Value> const& productItf) {
        return static_cast<Product<Mode, Value> const&>(productItf);
    };
}

#include "Product.impl.h"

#endif


