#include "Measurements.h"

namespace meas {
    
    std::int64_t reduce_steps(std::int64_t steps, All) {
        mpi::reduce<mpi::op::sum>(steps, mpi::master); return steps;
    };
    
    std::int64_t reduce_steps(std::int64_t steps, Jackknife) {
        auto temp = steps; mpi::all_reduce<mpi::op::sum>(temp); return temp - steps;
    };
    
    std::int64_t reduce_steps(std::int64_t steps, Variance) {
        mpi::all_reduce<mpi::op::sum>(steps); return steps;
    };
    
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
        
    
    void restart(jsx::value const& jIn, jsx::value & jMeas, double const eta, int const samples){
        
        if(jIn.is<jsx::array_t>()) {
            if(!(jMeas.is<jsx::array_t>() && jMeas.size() == jIn.size())) jMeas = jsx::array_t(jIn.size());
            int index = 0; for(auto const& jEntry : jIn.array()) restart(jEntry, jMeas[index++], eta, samples);
            
        } else if (jIn.is("io::cvec")){ //for some reason, jsx::value.is<io::rvec>() doesn't recognize the io::cvec?)
            auto jData = jIn("io::cvec");
            std::vector<ut::complex> data = jsx::at<io::cvec>(jData);
            for (auto& x : data) x*=samples*eta;
            jMeas << meas::fix(data, samples);
            
        } else if (jIn.is("io::rvec")) { //jIn.is<io::rvec>()
            auto jData = jIn("io::rvec");
            std::vector<double> data = jsx::at<io::rvec>(jData);
            for (auto& x : data) x*=samples*eta;
            jMeas << meas::fix(data, samples);
            
        } else if (jIn.is<jsx::object_t>()){
            for (auto const& jEntry : jIn.object()){
                if (jEntry.first != "steps" and jEntry.first != "eta"){ //not real measurements -- WangLandau deals with them
                    if (jEntry.first != "expansion histogram") //the only vecvar I've encountered -- TODO: an automatic way to do this
                        restart(jEntry.second, jMeas[jEntry.first], eta, samples);
                }
            }
        }
        
    }

    void restart(jsx::value & jParams, jsx::value & jMeasurements){
        
        if(jParams.is("restart") and jParams("restart").boolean()){
            
            mpi::cout << "Initializing measurements from previous run" << std::endl;
            auto const names = cfg::Worm::get_names();
            
            //The measurements stored during ctqmc are not yet normalized by the relative size of the configuration spaces
            //or the average sign, but the reported measurements in params.meas.json are
            //so, we must go back to the raw measurement
            //this is the 1st part of refined -> raw (see reduce, e.g.)
            auto const pName = cfg::partition::Worm::name();
            auto jSign = jParams("measurements")(pName)("sign")("io::rvec");
            std::vector<double> sign = jsx::at<io::rvec>(jSign);
            auto const pEta = jParams("measurements")(pName)("eta").real64();
            auto const pSteps = jParams("measurements")(pName)("steps").int64();
            auto const signxZp = sign.at(0)*pSteps/pEta;
            
            for (auto & jIn : jParams("measurements").object()){ //loop through configuration spaces
                
                //2nd part of the refined -> raw
                auto const Zw = jIn.second("steps").int64()/jIn.second("eta").real64();
                
                //weight previous measurements of by # of samples taken
                std::int64_t sweep = 0;
                if (jIn.first == cfg::partition::Worm::name()){
                    std::int64_t const sweepB = jParams(jIn.first).is("sweepB") ? jParams(jIn.first)("sweepB").int64() : 250;
                    std::int64_t const sweepA = jParams(jIn.first).is("sweepA") ? jParams(jIn.first)("sweepA").int64() : 50;
                    sweep = std::max(sweepB, sweepA);
                    //This will slow convergence for the smaller of sweepB/sweepB (by default, the A-type partition obs)
                    //but this is better than preventing convergence for the other type.
                    //Consider using sweepB only if susceptibility bulla or green bulla set
                    //TODO: It would be better if params.meas.json output the sweep / var / fix of each observable so that restart could use it
                } else {
                    sweep = jParams(jIn.first).is("sweep") ? jParams(jIn.first)("sweep").int64() : 50;
                }
                
                std::size_t samples = jIn.second("steps").int64() / (mpi::number_of_workers() * sweep);
                
                restart(jIn.second, jMeasurements[jIn.first], signxZp/Zw, samples);
                
            }
            
            jParams["measurements"]=jsx::empty_t(); //free memory -- old measurements
            
            mpi::cout << "Done initializing measurements from previous run" << std::endl;
        }
    }

    //Each worker is not guarenteed to have visited the same subspaces of the larger  worm space
    //This is not only concerning from a convergence standpoint, but it also will cause MPI to hang
    //When we try to reduce the results. So, we add a measurement of zero when there is such a miss
    template <typename Value>
    void check_missing_tensor_elements(jsx::value const& jParams, jsx::value& jMeasurements){
        
        if (mpi::number_of_workers() > 1)
            for(auto& space : jMeasurements.object())
                if (jParams.is(space.first) and space.first != cfg::partition::Worm::name()){
                    if (jParams(space.first)("basis").string() == "matsubara")
                        check_missing_tensor_elements<cvecfix,ut::complex>(space.first, space.second);
                    else if (jParams(space.first)("basis").string() == "legendre")
                        check_missing_tensor_elements<Vector<Value,Fix>,Value>(space.first, space.second);
                    
                }
        
    }
    
    template void check_missing_tensor_elements<double>(jsx::value const& jParams, jsx::value& jMeasurements);

    template void check_missing_tensor_elements<std::complex<double>>(jsx::value const& jParams, jsx::value& jMeasurements);
    
}

