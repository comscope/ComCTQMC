#ifndef INCLUDE_MEASUREMENTS_MEASUREMENTS_IMPL_H
#define INCLUDE_MEASUREMENTS_MEASUREMENTS_IMPL_H

#include "Measurements.h"

namespace meas {


    template<typename T, typename M>
    void Vector<T,M>::add(std::vector<T> const& val, std::int64_t samples) {
        resize_add(val, data_, M()); samples_ += samples; for(std::size_t n = 0; n < val.size(); ++n) data_[n] += val[n];
    };

    template<typename T, typename M>
    jsx::value Vector<T,M>::reduce(double fact, All, bool b64) const {
        auto samples = samples_;  mpi::barrier(); mpi::reduce<mpi::op::sum>(samples, mpi::master);
        auto data = data_;  resize_reduce(data, M());  mpi::barrier(); mpi::reduce<mpi::op::sum>(data, mpi::master);
        
        if(mpi::rank() == mpi::master) {
            if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
            for(auto& x : data) x *= fact/samples;
            
            data.b64() = b64;  return std::move(data);
        } else
            return jsx::null_t();
    };

    template<typename T, typename M>
    jsx::value Vector<T,M>::reduce(double fact, Jackknife, bool b64) const {
        auto samples = samples_; mpi::barrier(); mpi::all_reduce<mpi::op::sum>(samples);  samples -= samples_;
        
        if(!samples) throw std::runtime_error(name() + "::write: no measurements taken !");  //Scheisse das sštt eh nie passiere.
        
        auto data = data_;  resize_reduce(data, M()); mpi::barrier(); mpi::all_reduce<mpi::op::sum>(data);
        for(std::size_t i = 0; i < data_.size(); ++i) data[i] -= data_[i];
        for(auto& x : data) x *= fact/samples;
        
        data.b64() = b64; return std::move(data);
    };

    template<typename T, typename M>
    void Vector<T,M>::write(jsx::value& dest) const {
        throw std::runtime_error(name() + "::write: not implemented");
    };

    template<typename T, typename M>
    std::size_t Vector<T,M>::size() const {return data_.size();}
    template<typename T, typename M>
    T Vector<T,M>::at(int const i) const {return data_[i];}

    template<typename T, typename M>
    void Vector<T,M>::resize_add(io::Vector<T> const& val, io::Vector<T>& data, Fix) {
        if(data.size() == 0) data.resize(val.size(), .0);
        if(data.size() != val.size()) throw std::runtime_error(name() + "::add: missmatch in array size!");
    };
    template<typename T, typename M>
    void Vector<T,M>::resize_add(io::Vector<T> const& val, io::Vector<T>& data, Var) {
        if(val.size() > data.size()) data.resize(val.size(), .0);
    };

    template<typename T, typename M>
    void Vector<T,M>::resize_reduce(io::Vector<T>& data, Fix) {
    };
    template<typename T, typename M>
    void Vector<T,M>::resize_reduce(io::Vector<T>& data, Var) {
        auto size = data.size(); mpi::all_reduce<mpi::op::max>(size); data.resize(size, .0);
    };




    //--------------------------------------------------------------------------------------------------------------------------------


    template  <typename measure_type, typename Value>
    void check_missing_tensor_elements(std::string const name, jsx::value& jIn){
        
        for (auto& jMeas : jIn.object()){ //eta, steps, static, dynamic
            if (jMeas.second.is<jsx::object_t>()){ //skip "eta" and "steps"
                
                //First, create a string which is the concatenated list of tensory elements i_j_k_l
                //Also, get size of the measurement vector
                std::string concat_of_entries = "";
                std::vector<std::string> list_of_entries;
                std::size_t size=0;
                for (auto& jEntry : jMeas.second.object()){ //i_j_k_l's or i_j's
                    auto temp = jEntry.second.at<measure_type>().at(0); //if an entry is 0, then it will disapear during post processing...
                    if (std::abs(temp) > 1.e-14){
                        concat_of_entries+=jEntry.first+" ";
                        list_of_entries.push_back(jEntry.first);
                    }
                    size = jEntry.second.at<measure_type>().size(); //TODO: This should be io::Vector<Value> if in legendre basis...
                }
                
                //Pad the strings so that this list is the same length across all workers
                auto length_of_entries = concat_of_entries.size();
                mpi::all_reduce<mpi::op::max>(length_of_entries);
                concat_of_entries.resize(length_of_entries, ' ');
                
                //Also, set a vector which represents that no measurements have been taken.
                mpi::all_reduce<mpi::op::max>(size);
                std::vector<Value> const no_meas(size, 1.e-7); //TODO: This should be Value if in legendre basis...
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

}

#endif
