#include "Scheduler.h"

namespace mch {
    
    Scheduler::Scheduler(bool thermalised, Phase phase) : thermalised_(thermalised), phase_(phase) {};
    
    TimeScheduler::TimeScheduler(std::int64_t duration, bool thermalised, Phase phase, int64_t overtime) :
        Scheduler(thermalised, phase),
        duration_(60.*duration - overtime),
        start_(std::chrono::steady_clock::now()) { };
    
    std::int64_t TimeScheduler::overtime() const { return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() - duration_; }
    
        
    bool TimeScheduler::done() {
        return !(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() < duration_);
    };
    
    StepsScheduler::StepsScheduler(std::int64_t stop, bool thermalised, Phase phase) :
        Scheduler(thermalised, phase),
        stop_(stop),
        steps_(0) { };
    
    bool StepsScheduler::done() {
        return ++steps_ >= stop_;
    };
    
}











