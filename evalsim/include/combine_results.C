#include "combine_results.h"

namespace evalsim{

void MeasurementCombiner::combine(jsx::value const& jParams, std::string const name){
    auto dirs = jParams("directories").array();
    bool first = true;
    for (auto & d : dirs){
        
        auto dir_name = d.string();
        auto const jMeasurements = mpi::read(dir_name + name);
        
        if (first){ data_ = jMeasurements;}
        
        for (auto & jIn : jMeasurements.object()){
            add(jIn.second, data_[jIn.first], jIn.second("steps").real64(), first);
        }
        
        first = false;
    }
    
    for (auto & jIn : data_.object()){
        add(jIn.second, data_[jIn.first], 1.0/jIn.second("steps").real64(), true);
    }
}

void MeasurementCombiner::add(jsx::value const& jIn, jsx::value & jMeas, double const factor, bool const first){
    if(jIn.is<jsx::array_t>()) {
        if(!(jMeas.is<jsx::array_t>() && jMeas.size() == jIn.size())) jMeas = jsx::array_t(jIn.size());
        int index = 0; for(auto const& jEntry : jIn.array()) add(jEntry, jMeas[index++], factor, first);
        
    } else if (jIn.is("io::cvec")){ //for some reason, jsx::value.is<io::rvec>() doesn't recognize the io::cvec?)
        
        auto cvec1 = jIn("io::cvec");
        auto data1 = jsx::at<io::cvec>(cvec1);
        auto cvec2 = jMeas("io::cvec");
        auto data2 = jsx::at<io::cvec>(cvec2);
        
        for (int i=0; i<data1.size(); i++){
            if (first)
                data2[i] = data1[i]*factor;
            else
                data2[i] += data1[i]*factor;
        }
        jMeas["io::cvec"] = data2;
        
    } else if (jIn.is("io::rvec")) { //jIn.is<io::rvec>()
        
        auto rvec1 = jIn("io::rvec");
        auto data1 = jsx::at<io::rvec>(rvec1);
        auto rvec2 = jMeas("io::rvec");
        auto data2 = jsx::at<io::rvec>(rvec2);
        
        for (int i=0; i<data1.size(); i++){
            if (first)
                data2[i] = data1[i]*factor;
            else
                data2[i] += data1[i]*factor;
        }
        jMeas["io::rvec"] = data2;
        
    } else if (jIn.is<jsx::object_t>()){
        for (auto const& jEntry : jIn.object()){
            if (jEntry.first == "steps"){
                if (!first){
                    auto step1 = jMeas("steps").int64();
                    auto step2 = jIn("steps").int64();
                    jMeas["steps"] = step1 + step2;
                }
            }
            else if (jEntry.first != "eta"){ //not real measurements
                add(jEntry.second, jMeas[jEntry.first], factor, first);
            }
        }
    }
}
}
