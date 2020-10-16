#ifndef CTQMC_INCLUDE_MARKOVCHAIN_SCHEDULER_H
#define CTQMC_INCLUDE_MARKOVCHAIN_SCHEDULER_H

#include <ratio>
#include <chrono>
#include <ctime>
#include <tuple>
#include <random>


namespace mch {
    
    enum class Phase { Initialize, Step, Sample, Finalize };
    
    struct Scheduler {
        Scheduler() = delete;
        Scheduler(bool thermalised, Phase phase);
        bool thermalised() const { return thermalised_;};
        Phase& phase() { return phase_;}
        virtual std::int64_t overtime() const { return 0; }
        virtual bool done() { return true;};
        virtual ~Scheduler() = default;
    protected:
        bool const thermalised_;
        Phase phase_;
    };
    
    struct TimeScheduler : Scheduler {
        TimeScheduler(std::int64_t duration, bool thermalised, Phase phase, int64_t overtime = 0);
        ~TimeScheduler() = default;
        
        std::int64_t overtime() const;
        
        bool done();
        
    private:
        double const duration_;
        std::chrono::steady_clock::time_point const start_;
    };
    
    struct StepsScheduler : Scheduler {
        StepsScheduler(std::int64_t stop, bool thermalised, Phase phase);
        ~StepsScheduler() = default;
        
        bool done();
        
    private:
        std::int64_t const stop_;
        std::int64_t steps_;
    };
    
}

#include "Scheduler.impl.h"


#endif











