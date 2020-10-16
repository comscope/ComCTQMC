#ifndef INCLUDE_ATOMIC_GENERATE_IMPL_H
#define INCLUDE_ATOMIC_GENERATE_IMPL_H

#include "Generate.h"

namespace ga {
    

    
    template<Order order, typename Value>
    jsx::value get_observable(Tensor<Value> const& tensor,
                              BlockStates const& blockStates,
                              bool const throw_error) //if false, instead return jsx::empty on error
    {
        jsx::value jObservable = jsx::array_t(blockStates.size());

        State block_label = 0;
        for(auto const& blockState : blockStates) {
            jObservable[block_label]["target"] = jsx::int64_t(block_label);
            io::Matrix<Value> matrix(blockState.size(), blockState.size());
            
            State state_indexJ = 0;
            for(auto const & state : blockState) {
                FlavorState const stateJ(state);
                
                for(int f1 = 0; f1 < tensor.N(); ++f1)
                    for(int f2 = 0; f2 < tensor.N(); ++f2)
                        if(std::abs(tensor.t(f1, f2)) > 1.e-14) {
                            FlavorState const stateI = psiDagg(f1, psi(f2, stateJ));
                            
                            if(stateI.sign() != 0) {
                                auto it = std::find(blockState.begin(), blockState.end(), stateI.state());
                                
                                if(it == blockState.end()){
                                    if (throw_error) throw std::runtime_error("Something is wrong with the partitioning of the states");
                                    else {mpi::cout << " not a good observable; removing from list ... " ; return jsx::empty();}
                                }
                                
                                
                                matrix(it - blockState.begin(), state_indexJ) += tensor.t(f1, f2)*static_cast<double>(stateI.sign());
                            }
                        }
                
                for(int f1 = 0; f1 < tensor.N(); ++f1)
                    for(int f2 = 0; f2 < tensor.N(); ++f2)
                        for(int f3 = 0; f3 < tensor.N(); ++f3)
                            for(int f4 = 0; f4 < tensor.N(); ++f4)
                                if(std::abs(tensor.V(f1, f2, f3, f4)) > 1.e-14) {
                                    FlavorState const stateI = order == Order::normal ? psiDagg(f1, psiDagg(f2, psi(f3, psi(f4, stateJ)))) : psiDagg(f1, psi(f2, psiDagg(f3, psi(f4, stateJ))));
                                    
                                    if(stateI.sign() != 0) {
                                        auto it = std::find(blockState.begin(), blockState.end(), stateI.state());
                                        
                                        if(it == blockState.end()){
                                            if (throw_error) throw std::runtime_error("Something is wrong with the partitioning of the states");
                                            else {mpi::cout << " not a good observable; removing from list ... " ; return jsx::empty();}
                                        }
                                        
                                        matrix(it - blockState.begin(), state_indexJ) += tensor.V(f1, f2, f3, f4)*static_cast<double>(stateI.sign());
                                    }
                                }
                ++state_indexJ;
            }
            jObservable[block_label]["matrix"] = std::move(matrix);
        
            ++block_label;
        }
        
        return jObservable;
    };

    
};

#endif







