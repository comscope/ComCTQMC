#ifndef INCLUDE_MEASUREMENTS_MEASUREMENTS_H
#define INCLUDE_MEASUREMENTS_MEASUREMENTS_H

#include <vector>
#include <complex>
#include <fstream>
#include <valarray>
#include <cstring>
#include <random>

#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../io/Vector.h"
#include "../../ctqmc/include/config/Worms.h"

//Achtung: es kann sein dass gewisse observabeln nicht gespeichert wurden, c.f. MonteCarlo.h

namespace meas {
    
    struct All {}; struct Jackknife {};
    
    
    struct Fix {};  struct Var {};

    template<typename T, typename M>
    struct Vector {
        inline static std::string name() { return name(T(), M());};
        
        Vector() = default;
        Vector(Vector const&) = default;
        Vector(Vector&&) = default;
        Vector& operator=(Vector const&) = default;
        Vector& operator=(Vector&&) = default;
        ~Vector() = default;
        
        void add(std::vector<T> const& val, std::int64_t samples) {
            resize_add(val, data_, M()); samples_ += samples; for(std::size_t n = 0; n < val.size(); ++n) data_[n] += val[n];
        };
        
        jsx::value reduce(double fact, All, bool b64) const {
            auto samples = samples_;  mpi::reduce<mpi::op::sum>(samples, mpi::master);
            auto data = data_;  resize_reduce(data, M());  mpi::reduce<mpi::op::sum>(data, mpi::master);
            
            if(mpi::rank() == mpi::master) {
                if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
                for(auto& x : data) x *= fact/samples;
                data.b64() = b64;  return std::move(data);
            } else
                return jsx::null_t();
        };
        
        jsx::value reduce(double fact, Jackknife, bool b64) const {
            auto samples = samples_;  mpi::all_reduce<mpi::op::sum>(samples);  samples -= samples_;
            
            if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
            
            auto data = data_;  resize_reduce(data, M());  mpi::all_reduce<mpi::op::sum>(data);
            for(std::size_t i = 0; i < data_.size(); ++i) data[i] -= data_[i];
            for(auto& x : data) x *= fact/samples;
            
            data.b64() = b64; return std::move(data);
        };
        
        void write(jsx::value& dest) const {
            throw std::runtime_error(name() + "::write: not implemented");
        };
        
        std::size_t size() const {return data_.size();}
        T at(int const i) const {return data_[i];}
        
    private:
        
        std::int64_t samples_ = 0;
        io::Vector<T> data_;
        
        static std::string name(double const&, Fix)               { return "meas::rvecfix"; };
        static std::string name(std::complex<double> const&, Fix) { return "meas::cvecfix"; };
        
        static std::string name(double const&, Var)               { return "meas::rvecvar"; };
        static std::string name(std::complex<double> const&, Var) { return "meas::cvecvar"; };
        
        static void resize_add(io::Vector<T> const& val, io::Vector<T>& data, Fix) {
            if(data.size() == 0) data.resize(val.size(), .0);
            if(data.size() != val.size()) throw std::runtime_error(name() + "::add: missmatch in array size!");
        };
        static void resize_add(io::Vector<T> const& val, io::Vector<T>& data, Var) {
            if(val.size() > data.size()) data.resize(val.size(), .0);
        };
        
        static void resize_reduce(io::Vector<T>& data, Fix) {
        };
        static void resize_reduce(io::Vector<T>& data, Var) {
            auto size = data.size(); mpi::all_reduce<mpi::op::max>(size); data.resize(size, .0);
        };
    };

    
    template<typename T, typename M>
    struct Sample {
        T const& value; std::int64_t samples;
    };
    
    
    template<typename T> inline Sample<T, Fix> fix(T const& value, std::int64_t samples) { return {value, samples};};
    template<typename T> inline Sample<T, Var> var(T const& value, std::int64_t samples) { return {value, samples};};
    
    
    template<typename T, typename M>
    inline void operator<<(jsx::value& lhs, Sample<T, M> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = Vector<T, M>();
        lhs.at<Vector<T, M>>().add(std::vector<T>(1, rhs.value), rhs.samples);
    }
    
    template<typename T, typename M>
    inline void operator<<(jsx::value& lhs, Sample<std::vector<T>, M> const& rhs) {
        if(lhs.is<jsx::empty_t>()) lhs = Vector<T, M>();
        lhs.at<Vector<T, M>>().add(rhs.value, rhs.samples);
    }

    
    using rvecfix = Vector<double, Fix>;  using cvecfix = Vector<std::complex<double>, Fix>;
    using rvecvar = Vector<double, Var>;  using cvecvar = Vector<std::complex<double>, Var>;
    
    
    //--------------------------------------------------------------------------------------------------------------------------------
    
    
    std::int64_t reduce_steps(std::int64_t steps, All) {
        mpi::reduce<mpi::op::sum>(steps, mpi::master); return steps;
    };
    
    std::int64_t reduce_steps(std::int64_t steps, Jackknife) {
        auto temp = steps; mpi::all_reduce<mpi::op::sum>(temp); return temp - steps;
    };
    
    
    template<typename E>
    inline void reduce(jsx::value& jOut, double fact, jsx::value const& jIn, E, bool b64) {
        if(jIn.is<rvecfix>())
            jOut = jIn.at<rvecfix>().reduce(fact, E(), b64);
        else if(jIn.is<cvecfix>())
            jOut = jIn.at<cvecfix>().reduce(fact, E(), b64);
        else if(jIn.is<rvecvar>())
            jOut = jIn.at<rvecvar>().reduce(fact, E(), b64);
        else if(jIn.is<cvecvar>())
            jOut = jIn.at<cvecvar>().reduce(fact, E(), b64);
        else if(jIn.is<jsx::object_t>()) {
            for(auto& jEntry : jIn.object()) reduce(jOut[jEntry.first], fact, jEntry.second, E(), b64);
        } else if(jIn.is<jsx::array_t>()) {
            if(!(jOut.is<jsx::array_t>() && jOut.size() == jIn.size())) jOut = jsx::array_t(jIn.size());
            int index = 0; for(auto& jEntry : jIn.array()) reduce(jOut[index++], fact, jEntry, E(), b64);
        } else
            jOut = jIn;
    }
    
    
    template<typename E>
    inline void reduce(jsx::value& jOut, jsx::value const& jIn, jsx::value const& jEtas, E, bool b64) {
        auto const pName = cfg::partition::Worm::name();
        
        jsx::value const jSign = jIn(pName)("sign").at<rvecfix>().reduce(1., E(), b64);
        auto const pSteps = reduce_steps(jIn(pName)("steps").int64(), E());
        auto const signxZp = (jSign.is<jsx::null_t>() ? 1. : jsx::at<io::rvec>(jSign).at(0))*pSteps/jEtas(pName).real64();
        
        if(!jOut.is<jsx::object_t>()) jOut = jsx::object_t();
        
        for(auto& jWorm : jIn.object()) {
            auto const wSteps = reduce_steps(jWorm.second("steps").int64(), E());
            auto const Zw = wSteps/jEtas(jWorm.first).real64();
            
            reduce(jOut[jWorm.first], Zw/signxZp, jWorm.second, E(), b64);
            
            jOut[jWorm.first]["steps"] = wSteps;
        }
        
        jOut[pName]["sign"] = jSign;
    }

std::vector<std::string> split_by_char(std::string const& string, char const c){
    std::stringstream temp(string);
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(temp, segment, c))
    {
       seglist.push_back(segment);
    }
    
    return seglist;
}

    void check_missing_tensor_elements(std::string const name, jsx::value& jIn){
        
        for (auto& jMeas : jIn.object()){ //eta, steps, static, dynamic
            if (jMeas.second.is<jsx::object_t>()){ //skip "eta" and "steps"
                
                //First, create a string which is the concatenated list of tensory elements i_j_k_l
                //Also, get size of the measurement vector
                std::string concat_of_entries = "";
                std::vector<std::string> list_of_entries;
                std::size_t size=0;
                for (auto& jEntry : jMeas.second.object()){ //i_j_k_l's or i_j's
                    auto temp = jEntry.second.at<cvecfix>().at(0); //if an entry is 0, then it will disapear during post processing...
                    if (std::abs(temp) > 1.e-14){
                        concat_of_entries+=jEntry.first+" ";
                        list_of_entries.push_back(jEntry.first);
                    }
                    size = jEntry.second.at<cvecfix>().size(); //TODO: This should be io::Vector<Value> if in legendre basis...
                }
                
                //Pad the strings so that this list is the same length across all workers
                auto length_of_entries = concat_of_entries.size();
                mpi::all_reduce<mpi::op::max>(length_of_entries);
                concat_of_entries.resize(length_of_entries, ' ');
                
                //Also, set a vector which represents that no measurements have been taken.
                mpi::all_reduce<mpi::op::max>(size);
                std::vector<ut::complex> const no_meas(size, 1.e-7); //TODO: This should be Value if in legendre basis...
                std::int64_t const no_samples = 1; //Causes errors if it's zero...
                
                //Concatenate list of entries  from all workers (on master rank)
                mpi::gather(concat_of_entries, mpi::master);
                
                //TODO: not actually gathering unique elements
                //Gather only the unique elements from this list
                std::string unique_concat_of_entries="";
                if (mpi::rank() == mpi::master){
                    concat_of_entries.pop_back();
                    
                    auto unique_list_of_entries = split_by_char(concat_of_entries,' ');
                    std::sort(unique_list_of_entries.begin(),unique_list_of_entries.end());
                    
                    std::vector<std::string>::iterator it;
                    it = std::unique(unique_list_of_entries.begin(), unique_list_of_entries.end());
                    
                    auto const number_of_unique_entries = std::distance(unique_list_of_entries.begin(),it);
                    unique_list_of_entries.resize(number_of_unique_entries);
                    
                    for (auto const& entry : unique_list_of_entries){
                        unique_concat_of_entries += entry + " ";
                    }
                    unique_concat_of_entries.pop_back();
                }
                
                
                //Broadcast this list back to all workers
                size = unique_concat_of_entries.size();
                mpi::bcast(size, mpi::master);
                unique_concat_of_entries.resize(size);
                mpi::bcast(unique_concat_of_entries, mpi::master);
                
                auto const unique_list_of_entries = split_by_char(unique_concat_of_entries,' ');
                
                //Check that this list is the same as the original list (i.e., each worker has at least had one step in each space all other workers have  sampled)
                int missing = 0;
                for (auto  const& entry : unique_list_of_entries){
                    if (entry != "" and entry != " "){ //TODO: get rid of need for this check
                        auto const check = std::find(list_of_entries.begin(),list_of_entries.end(),entry);
                        
                        //if an entry is not in the workers list of visited spaces, initialize this space
                        //And keep track of how many spaces are not hit
                        if (check == list_of_entries.end()){
                            missing++;
                            jMeas.second[entry] << fix(no_meas, no_samples);
                        }
                    }
                }
            
                if (missing)
                    mpi::cout << "Warning: Each worker did not sample all sub-spaces of worm space " << name << ": " << missing << " misses" << std::endl;
                
            }
        }
    }
    

    //Each worker is not guarenteed to have visited the same subspaces of the larger  worm space
    //This is not only concerning from a convergence standpoint, but it also will cause MPI to hang
    //When we try to reduce the results. So, we add a measurement of zero when there is such a miss
    void check_missing_tensor_elements(jsx::value const& jParams, jsx::value& jMeasurements){

        if (mpi::number_of_workers() > 1)
            for(auto& space : jMeasurements.object())
                if (jParams.is(space.first) and space.first != cfg::partition::Worm::name()){
                    check_missing_tensor_elements(space.first, space.second);
                }
        
    }

    void restart(jsx::value const& jParams, jsx::value & jInput, jsx::value & jMeasurements){
        if(1){
            
            mpi::cout << "Initializing measurements from previous run" << std::endl;
            auto const names = cfg::Worm::get_names();
            
            for (auto & meas : jInput.object()){
                
                auto samples = jInput[meas.first]["steps"].int64() / mpi::number_of_workers();
                if (mpi::rank() == mpi::master) samples += jInput[meas.first]["steps"].int64() % mpi::number_of_workers();
                
                //TODO: test
                //TODO: Add partition space restart
                auto it = std::find(names.begin(), names.end(), meas.first);
                if(it != names.end() and meas.first == cfg::partition::Worm::name()){
                    
                    for (auto & type : meas.second.object()){ //static vs dynamic
                        for (auto & func : type.second.object()){ //ij or ijkl of space
                                
                            //need to do this for the other types
                            if (func.second.is<io::cvec>()){
                                std::vector<ut::complex> data = jsx::at<io::cvec>(func.second);
                                for (auto& x : data) x*=samples;
                                jMeasurements(meas.first)(type.first)(func.first) << meas::fix(data, samples);
                                
                            } else if (func.second.is<io::rvec>()) {
                                std::vector<double> data =jsx::at<io::rvec>(func.second);
                                for (auto& x : data) x*=samples;
                                jMeasurements(meas.first)(type.first)(func.first) << meas::fix(data, samples);
                            }
                        
                        }
                    }
                }
                
            }
         
            mpi::cout << "Done initializing measurements from previous run" << std::endl;
        }
    }

}




#endif
