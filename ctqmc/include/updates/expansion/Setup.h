#ifndef CTQMC_INCLUDE_UPDATES_EXPANSION_SETUP_H
#define CTQMC_INCLUDE_UPDATES_EXPANSION_SETUP_H

#include <vector>
#include <array>

#include "InsertErasePairCsc.h"
#include "InsertErasePairStd.h"
#include "InsertEraseQuadStd.h"
#include "GlobalSpinFlip.h"
#include "OpShift.h"

#include "../../markovchain/Update.h"
#include "../../markovchain/MarkovChain.h"

// Scheisse kack 4 operator updates immer no nit implementiert Himmelhergottsakramentzifixhallelujaleckstmiamarscheisseglumpvereckts !

namespace upd {
    
    namespace expansion {
        
        template<typename Space, typename Mode, typename Value>
        void setup_updates(jsx::value const& jParams, data::Data<Value> const& data, state::State<Value> const& state, mch::MarkovChain<Value>& markovChain, int const stream) {
            
            if (!stream)
                mpi::cout << "Setting expansion updates for " + Space::name() + " space ... ";
            
            markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< InsertCsc<Space>, Mode, Value >(1., jParams, data)),
                            mch::unique_update_ptr<Value>(new upd::Generic< EraseCsc<Space>, Mode, Value >(1., jParams, data)));
            
            if (jParams("quad insert").boolean()){
                        markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< QuadInsertStd<Space>, Mode, Value >(1., jParams, data)),
                                        mch::unique_update_ptr<Value>(new upd::Generic< QuadEraseStd<Space>, Mode, Value >(1., jParams, data)));
            }
            
            if (jParams("spin flip").real64() > 0){
                markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< SpinFlip<Space>, Mode, Value >(jParams("spin flip").real64(), jParams, data)),
                                mch::unique_update_ptr<Value>(new upd::Generic< SpinFlip<Space>, Mode, Value >(jParams("spin flip").real64(), jParams, data)));
            }
            
            if (jParams("shift").boolean()){
                markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< RShift<Space>, Mode, Value >(1.0, jParams, data)),
                                mch::unique_update_ptr<Value>(new upd::Generic< RShift<Space>, Mode, Value >(1.0, jParams, data)));
                
                markovChain.add(mch::unique_update_ptr<Value>(new upd::Generic< LShift<Space>, Mode, Value >(1.0, jParams, data)),
                                mch::unique_update_ptr<Value>(new upd::Generic< LShift<Space>, Mode, Value >(1.0, jParams, data)));
            }
            
            if (!stream)
                mpi::cout << "Ok" << std::endl;
            
        };
        
    }
    
}

#endif
