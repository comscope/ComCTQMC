#include "Initialize.h"

namespace params {
/*
This is an effort to more sensibly fill out and check the parameter file supplied to ComCTQMC
The goal is to eliminate all checks (jParams.is() ? jParams() : default)
And provide more sensible errors when users mess up the input
It would be nice to do this in a more automatic fashion...
*/

    void check_type_against_default(std::string const& entry, jsx::value const& val, jsx::value const& default_val){
        
        if (default_val.is<jsx::boolean_t>() and !val.is<jsx::boolean_t>())
            throw std::runtime_error("Expected boolean for parameter named: " + entry);
        
        //allow promotion of integers to real
        if (default_val.is<jsx::real64_t>() and (!val.is<jsx::real64_t>() and !val.is<jsx::int64_t>()))
            throw std::runtime_error("Expected real number for parameter named: " + entry);
        
        //allow demotion of reals to integers -- can occasionally cause issues although jsonx tries to cast int to real
        if (default_val.is<jsx::int64_t>() and (!val.is<jsx::real64_t>() and !val.is<jsx::int64_t>()))
            throw std::runtime_error("Expected integer for parameter named: " + entry);
        
        if (default_val.is<jsx::string_t>() and !val.is<jsx::string_t>())
            throw std::runtime_error("Expected string for parameter named: " + entry);
        
        if (default_val.is<jsx::array_t>() and !val.is<jsx::array_t>()){
            throw std::runtime_error("Expected array for parameter named: " + entry);
            if (default_val.size() == val.size())
                for (auto const& x : val.array())
                    check_type_against_default(entry,x,default_val.array()[0]);
            }
        
        if (default_val.is<jsx::object_t>() and !val.is<jsx::object_t>()){
            throw std::runtime_error("Expected array for parameter named: " + entry);
            if (default_val.size() == val.size()){
                for (auto const& x : default_val.object())
                    check_type_against_default(entry + " -> " + x.first, val(x.first), x.second);
            }
        }
    }

    void check_required_input(jsx::value const& jParams){
         
        MainRequired required;
        for (auto const& entry : required.get().object()){
            
            if (!jParams.is(entry.first))
                throw std::runtime_error("Parameter file is missing the key " + entry.first);
            
            if (!entry.second.is<jsx::empty_t>())
                for (auto const& sub_entry : entry.second.array())
                    if (!jParams(entry.first).is(sub_entry.string()))
                        throw std::runtime_error("Parameter file is missing the key " + sub_entry.string() + " in block " + entry.first);
            
        }
        
        //check data-types of all required imput except the two-body part which doesn't fit this paradigm
        required.init_types();
        for (auto const& entry : required.get().object()){
            check_type_against_default(entry.first, jParams(entry.first), entry.second);
        }
        
        //check data types of two body input -- basis gets checked in opt::basis?
        jsx::value const& jTwoBody = jParams("hloc")("two body");
        if (jTwoBody.is<jsx::object_t>() and !jTwoBody.is("imag")){
            
            if (!jParams.is("basis"))
                throw std::runtime_error("If the two body tensor is not explicity defined, you must supply the one-particle basis");
            
        }
             
         
     }


    void add_default_worms( jsx::value & jParams, AllDefaults const& defaults) {
        
        for (auto const& worm : defaults.listOfDefaultWorms){
            jParams[worm] = jsx::object_t();
            jParams[worm]["meas"] = jsx::array_t({"imprsum"});
        }
        
    }
 
    void set_default_values(jsx::value & jParams){
        
        AllDefaults defaults(jParams);
        
        for (auto const& entry : defaults.mainDefaults.get().object()){
            if (!jParams.is(entry.first)) jParams[entry.first] = entry.second;
            else check_type_against_default("Main block -> " + entry.first, jParams[entry.first], entry.second);
        }
        
        for (auto const& entry : defaults.partitionSpaceDefaults.get().object()){
            if (!jParams(cfg::partition::Worm::name()).is(entry.first)) jParams[cfg::partition::Worm::name()][entry.first] = entry.second;
            else check_type_against_default(cfg::partition::Worm::name() + " -> " + entry.first, jParams[cfg::partition::Worm::name()][entry.first], entry.second);
        }
        
        cfg::for_each_type<cfg::Worm>::apply(worm_defaults_functor(), jParams, defaults);
    }
         
    bool is_complex_in_cvec(io::cvec const& vec){
        
        bool check = false;
        
        for (auto const& x : vec)
            check = check or std::abs(x.imag()) > 1e-14;
        
        return check;
    }

    bool is_complex_in_cmat(io::prettycmat const& mat){
        
        bool check = false;
        
        for (int i=0; i<mat.I(); i++)
            for (int j=0; j<mat.J(); j++)
                check = check or std::abs(mat(i,j).imag()) > 1e-14;
        
        return check;
    }

    bool validate_complex(jsx::value const& jParams){
        
        if (jParams("complex").boolean()){
            
            bool is_complex = false;
            
            if (jParams.is("basis") and jParams("basis").is("transformation")){
                
                auto jTransformation = jParams("basis")("transformation");
                auto const cmat = jsx::at<io::prettycmat>( jTransformation );
                is_complex = is_complex or is_complex_in_cmat( cmat );
                
            }
            
            if (jParams("hloc").is("two body") and jParams("hloc")("two body").is<io::cvec>()){
                
                auto jTwoBody = jParams("hloc")("two body");
                auto const cvec = jsx::at<io::cvec>( jTwoBody );
                is_complex = is_complex or is_complex_in_cvec( cvec );
                
            }

            auto jOneBody = jParams("hloc")("one body");
            auto const cmat = jsx::at<io::prettycmat>( jOneBody );
            is_complex = is_complex or is_complex_in_cmat( cmat );
            
            return is_complex;
            
        } else {
            return true;
        }
    }

    io::rvec cvec_to_rvec(io::cvec const& cv){
        io::rvec rv(cv.size());
        
        for (int i=0; i<cv.size(); i++) rv[i] = cv[i].real();
        
        return rv;
    }

    io::prettyrmat cmat_to_rmat(io::prettycmat const& cmat){
        io::prettyrmat rmat(cmat.I(),cmat.J());
        
        for (int i=0; i<cmat.I(); i++)
            for (int j=0; j<cmat.J(); j++)
                rmat(i,j) = cmat(i,j).real();
        
        return rmat;
    }

    void change_fake_complex_to_real(jsx::value & jParams){
        
        if (!validate_complex(jParams)){
            jParams["complex"] = false;
            
            if (jParams.is("basis") and jParams("basis").is("transformation")){
                
                auto jTransformation = jParams("basis")("transformation");
                auto const cmat = jsx::at<io::prettycmat>( jTransformation );
                jParams["basis"]["transformation"] = cmat_to_rmat(cmat);
                
            }
            
            if (jParams("hloc").is("two body") and jParams("hloc")("two body").is<io::cvec>()){
                
                auto jTwoBody = jParams("hloc")("two body");
                auto const cvec = jsx::at<io::cvec>( jTwoBody );
                jParams["hloc"]["two body"] = cvec_to_rvec(cvec);
                
            }
            
            auto jOneBody = jParams("hloc")("one body");
            auto const cmat = jsx::at<io::prettycmat>( jOneBody );
            jParams["hloc"]["one body"] = cmat_to_rmat(cmat);
            
        }
        
    }

    void initialize(jsx::value & jParams){
        
        change_fake_complex_to_real(jParams);
        check_required_input(jParams);
        set_default_values(jParams);
        mpi::write(jParams, "defaults.json");
        
    }


}

