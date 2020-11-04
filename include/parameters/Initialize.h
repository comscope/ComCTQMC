#ifndef INCLUDE_PARAMETERS_INITIALIZE_H
#define INCLUDE_PARAMETERS_INITIALIZE_H

#include <complex>

#include "Defaults.h"
#include "../JsonX.h"
#include "../mpi/Utilities.h"
#include "../../ctqmc/include/config/Worms.h"

namespace params {
/*
This is an effort to more sensibly fill out and check the parameter file supplied to ComCTQMC
The goal is to eliminate all checks (jParams.is() ? jParams() : default)
And provide more sensible errors when users mess up the input
It would be nice to do this in a more automatic fashion...
*/

//This is the main interface -- it will check the input jParams and fill it out with defaults
    void initialize(jsx::value & jParams);

    void check_type_against_default(std::string const& entry, jsx::value const& val, jsx::value const& default_val);

    void check_required_input(jsx::value const& jParams);

    struct worm_defaults_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value & jParams, AllDefaults const& defaults) const;
        
        void add_default_worms( jsx::value & jParams, AllDefaults const& defaults);
    };

    struct worm_cutoffs_functor {
        template<typename W>
        void operator()(ut::wrap<W> w, jsx::value & jParams, int const fermion_cutoff, int const boson_cutoff) const;
    };

    void set_default_values(jsx::value & jParams);

    bool is_complex_in_cvec(io::cvec const& vec);

    bool is_complex_in_cmat(io::prettycmat const& mat);

    bool validate_complex(jsx::value const& jParams);

    io::rvec cvec_to_rvec(io::cvec const& cv);

    io::prettyrmat cmat_to_rmat(io::prettycmat const& cmat);

    void change_fake_complex_to_real(jsx::value & jParams);
        
}

#include "Initialize.impl.h"

#endif
