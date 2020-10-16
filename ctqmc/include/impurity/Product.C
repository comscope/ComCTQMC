
#include "Product.h"
#include "../Algebra.h"


namespace imp {
    
    template<typename Mode, typename Value>
    Product<Mode,Value>::Product(jsx::value const& jParams, itf::EigenValues const& eig, itf::Operator<Value> const& ide, itf::Operators<Value> const& ops) :
    eig_(get<Mode>(eig)), ide_(get<Mode>(ide)), ops_(get<Mode>(ops)),
    urng_(std::mt19937(234), std::uniform_real_distribution<double>(.0, 1.)),
    prob_(jParams.is("skip-list probability") ? jParams("skip-list probability").real64() : .5),
    baseProb_(std::pow(prob_, (jParams.is("skip-list shift") ? jParams("skip-list shift").int64() : 0) + 1)),
    maxHeight_(50),
    height_(1), heightBackup_(1),
    size_(0), sizeBackup_(0),
    sign_(1),
    first_(new NodeType(         0, maxHeight_ + 1,   &ide_, -1, eig_)),
    last_( new NodeType(ut::KeyMax,              0, nullptr, -1, eig_)) {
        for(int l = 0; l <= height_; ++l) { first_->next[l] = last_; first_->entries[l] = size_ + 1;}
        first_->accept();
        
        all_.begin = new SectorNorm*[(maxHeight_ + 2)*eig.sectorNumber()];
        
        int maxDim = 0; for(int s = eig.sectorNumber(); s; --s) maxDim = std::max(eig_.at(s).dim(), maxDim);
        bufferA_.reset(new Matrix<Mode, Value>(maxDim*maxDim));
        bufferB_.reset(new Matrix<Mode, Value>(maxDim*maxDim));
        bufferC_.reset(new Matrix<Mode, Value>(maxDim*maxDim));
    };
    
    template<typename Mode, typename Value>
    Product<Mode,Value>::~Product() {
        reject();
        
        delete[] all_.begin;

        auto it = first_;
        while(it != last_) {
            auto temp = it;
            it = it->next[0];
            delete temp;
        }
        delete last_;
    };
    
    template<typename Mode, typename Value>
    bool Product<Mode,Value>::insert(ut::KeyType const key, int flavor) {
        return last() != insert_impl(key, random_height(), &ops_.at(flavor), flavor, eig_);
    };
    
    template<typename Mode, typename Value>
    typename Product<Mode, Value>::CAccessType Product<Mode,Value>::insert(ut::KeyType const key, itf::Operator<Value> const* op, int flavor) {
        return insert_impl(key, random_height(), &get<Mode>(*op), flavor, eig_);
    };
    
    template<typename Mode, typename Value>
    typename Product<Mode, Value>::CAccessType Product<Mode,Value>::insert(ut::KeyType const key, itf::Operator<Value> const* op, int flavor, int height) {
        if(height > maxHeight_ + 1) throw std::runtime_error("imp::Product::insert: invalid height");
        return insert_impl(key, height, &get<Mode>(*op), flavor, eig_);
    };
    
    template<typename Mode, typename Value>
    void Product<Mode,Value>::erase(ut::KeyType key) { erase_impl(key); };
    
    template<typename Mode, typename Value>
    std::vector<int> Product<Mode,Value>::map(CAccessType cbegin, int level, std::vector<int> const& sectors) {
        auto begin = AccessType(cbegin.node_); std::vector<int> result;

        begin->op(level)->assign(sectors, all_);
        begin.next(level) != begin.next(0) ? map_impl(begin, level, all_) : begin->prop()->add(all_);
        for(auto sec : sectors)
            result.push_back(cbegin->op(level)->map(sec).sector);

        return result;
    };
    
    template<typename Mode, typename Value>
    void Product<Mode,Value>::multiply(CAccessType cbegin, int level, int sec, itf::Batcher<Value>& batcher) {
        auto begin = AccessType(cbegin.node_); if(begin->op(level)->isMat(sec)) return;

        if(begin.next(level) != begin.next(0)) {
            if(size_ + 2 > ptrMat_.size()) { ptrMat_.resize(size_ + 2); ptrProp_.resize(size_ + 2);}
            multiply_impl(begin, level, sec,  ptrMat_.data(), ptrProp_.data(), bufferA_.get(), bufferB_.get(), batcher);
        } else {
            copyEvolveL(begin->op(level)->mat(sec, ide_.mat(sec).I()*ide_.mat(sec).J()), begin->prop()->at(sec), begin->op0->mat(sec), batcher);
            auto& map = begin->op(level)->set_map(sec); map.sector = begin->op0->map(sec).sector; norm(&map.norm, begin->op(level)->mat(sec), batcher);
        }
    };
    
    template<typename Mode, typename Value>
    int Product<Mode,Value>::accept() {
        heightBackup_ = height_;
        
        for(auto ptr : touched_) ptr->accept();
        for(auto ptr : inserted_) ptr->accept();
        for(auto ptr : erased_) delete ptr;
        
        touched_.clear(); inserted_.clear(); erased_.clear();

        sizeBackup_ = size_;
        
        int temp = sign_; sign_ = 1; return temp;
    };
    
    template<typename Mode, typename Value>
    void Product<Mode,Value>::reject() {
        height_ = heightBackup_;
        
        for(auto ptr : touched_) ptr->reject();
        for(auto ptr : inserted_) delete ptr;
        
        touched_.clear(); inserted_.clear(); erased_.clear();

        size_ = sizeBackup_;
        
        sign_ = 1;
    };
    
    template<typename Mode, typename Value>
    int Product<Mode,Value>::random_height() {
        int h = 1; double u = urng_();
        for(double p = baseProb_; u < p && h < maxHeight_; ++h, p *= prob_);
        return h;
    };
    
    template<typename Mode, typename Value>
    void Product<Mode,Value>::erase_impl(ut::KeyType const key) {
        if(!(0 < key && key < ut::KeyMax)) throw std::runtime_error("imp::Product::erase_impl: invalid key");
        
        std::vector<NodeType*> node(maxHeight_ + 1); NodeType* n = first_;
        int p = 0;
        
        for(int l = std::min(height_, maxHeight_); l >= 0; --l) {
            for(; key > n->next[l]->key; n = n->next[l]) p += n->entries[l];
            node[l] = n; --n->entries[l];
        }
        n = n->next[0];
        
        if(key != n->key) throw std::runtime_error("imp::Product::erase_impl: key not found !");
        
        erased_.push_back(n);
        ++p; --size_;
        
        for(int l = 0; l < n->height; ++l) {
            node[l]->next[l] = n->next[l];
            node[l]->entries[l] += n->entries[l];
        };
        
        if(node[0]->touch(0)) touched_.push_back(node[0]);
        for(int l = 1; l <= std::min(height_, maxHeight_); ++l)
            if(node[l - 1] != node[l])
                if(node[l]->touch(l)) touched_.push_back(node[l]);

        for(; height_ > 1 && first_->next[height_ - 1] == last_; --height_)
            ;
        
        sign_ *= p%2 ? -1 : 1;  //modulo hier ziemlich sicher nicht notwending ... guck standard wegen wraparound .... pass auf in sign() falls nicht notwending ...
    };
    
    template<typename Mode, typename Value>
    void Product<Mode,Value>::map_impl(AccessType it, int level, SectorNormPtrs& requested) {
        auto const last = it.next(level);
        while(it.next(level) == last) --level;

        do {
            if(it.next(level) != it.next(0)) {
                SectorNormPtrs missing;
                if(it->op(level)->missing(missing, requested)) map_impl(it, level, missing);
                it->op(level)->map(requested);
            } else {
                it->op0->map(requested);
                it->prop()->add(requested);
            }
            it = it.next(level);
        } while(it != last);
    }
    
    template<typename Mode, typename Value>
    void Product<Mode,Value>::multiply_impl(AccessType const begin, int const level, int const sec, Matrix<Mode, Value> const** const mat, Vector<Mode> const** const prop, Matrix<Mode, Value>* A, Matrix<Mode, Value>* B, itf::Batcher<Value>& batcher) {
        auto const end = begin.next(level);
        auto l = level; while(begin.next(l) == end) --l;

        auto s = sec; auto m = mat; auto p = prop;
        for(auto it = begin; it != end; it = it.next(l))
            if(it.next(l) != it.next(0)) {
                if(!it->op(l)->isMat(s)) multiply_impl(it, l, s, m, p, A, B, batcher);
                *m++ = &it->op(l)->mat(s); s = it->op(l)->map(s).sector; *p++ = nullptr;
            } else {
                *m++ = &it->op0->mat(s); s = it->op0->map(s).sector; *p++ = &it->prop()->at(s);
            }

        *m = nullptr; m = mat; p = prop;
        if(*p != nullptr) {
            copyEvolveL(*A, **p, **m, batcher);
            *m = A; std::swap(A, B);
        }
        while(*(++m + 1) != nullptr) {
            mult(*A, **m, **(m - 1), batcher);
            if(*++p != nullptr) evolveL(**p, *A, batcher);
            *m = A; std::swap(A, B);
        }
        mult(begin->op(level)->mat(sec, (*m)->I()*(*(m - 1))->J()), **m, **(m - 1), batcher);
        if(*++p != nullptr) evolveL(**p, begin->op(level)->mat(sec), batcher);
        
        auto& map = begin->op(level)->set_map(sec);
        map.sector = s; norm(&map.norm, begin->op(level)->mat(sec), batcher);
    };

    template struct Product<imp::Host, double>;
    template struct Product<imp::Host, ut::complex>;
    
#ifdef MAKE_GPU_ENABLED
    template struct Product<imp::Host, double>;
    template struct Product<imp::Host, ut::complex>;
#endif
    
}

