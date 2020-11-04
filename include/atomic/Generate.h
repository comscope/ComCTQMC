#ifndef INCLUDE_ATOMIC_GENERATE_H
#define INCLUDE_ATOMIC_GENERATE_H

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <set>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>

#include "Tensor.h"
#include "../JsonX.h"
#include "../io/Vector.h"
#include "../io/Matrix.h"
#include "../linalg/LinAlg.h"
#include "../mpi/Utilities.h"

namespace ga {
    //kack isch veralted ....
    
    struct FlavorState : public std::bitset<8*sizeof(unsigned long long)> {
        typedef unsigned long long State;
        
        FlavorState();
        explicit FlavorState(State const& state) : std::bitset<8*sizeof(State)>(state), sign_(1) {};
        FlavorState(FlavorState const& state);
        
        State state() const;
        int sign() const;
        int& sign();
    private:
        int sign_;
    };
    
    
    FlavorState psi(int flavor, FlavorState const& state);
        
    FlavorState psiDagg(int flavor, FlavorState const& state);
    
    typedef FlavorState::State State;
    typedef std::vector<State> States;
    typedef std::vector<States> BlockStates;
    
    struct Join {
        Join(std::size_t size);
        void join(State state1, State state2);
        std::size_t clean();
        State label(State state);
    private:
        std::vector<State> labels_;
        
        State find_representative(State label);
    };
    
    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    
    template<typename Value>
    void find_blocks(Tensor<Value> const& hloc, BlockStates& blockStates);

    enum class Order { alternating, normal };
    
    template<Order order, typename Value>
    jsx::value get_observable(Tensor<Value> const& tensor,
                              BlockStates const& blockStates,
                              bool const throw_error = true); //if false, instead return jsx::empty on error


    
    template<typename Value>
    jsx::value diagonalise(jsx::value& jHamiltonian);
    
    
    template<typename Value>
    void transform(jsx::value const& jTransformation,
                   jsx::value& jOperator);
    
    
    template<typename Value>
    jsx::value get_annihilation_operators(int N,
                                          BlockStates const& blockStates);

    io::rvec get_sector_qn(BlockStates const& blockStates,
                             std::vector<double> const& qn,
                           bool const throw_error = true);
    
    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    
    template<typename Value>
    jsx::value construct_hloc(jsx::value jTensors, bool b64 = true);
    
    
    //-----------------------------------------------------------------------------------------------------
    
    template<typename Value>
    jsx::value read_hloc(std::string const name);
    
    
    BlockStates get_block_states(jsx::value const& jHloc);
    
    
    //-----------------------------------------------------------------------------------------------------
    
    io::rvec construct_sector_qn(jsx::value const& jHloc, jsx::value jqn, bool throw_error = true);
    
    template<typename Value>
    jsx::value construct_annihilation_operators(jsx::value const& jHloc);

    template<typename Value>
    jsx::value construct_observable(jsx::value const& jHloc, jsx::value const& jTensors, bool const throw_error = true);
    
    template<typename Value>
    jsx::value construct_occupation_states(jsx::value const& jHloc);

};

#include "Generate.impl.h"

#endif







