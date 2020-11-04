
#include "ctqmc/include/Params.h"

#include "include/parameters/Initialize.h"

ut::Beta ut::beta;

template <typename T>
typename std::enable_if<!std::is_integral<T>::value, std::vector<T>>::type linspace(std::size_t const size, T const begin, T const end, bool endpoint = true){

    std::vector<T> v(size);
    
    if (size){
        
        auto const delta = (end - begin) / static_cast<double>(size - endpoint);
        
        v[0] = begin;
        for (std::size_t i = 1; i<size-endpoint; i++)
            v[i] = v[i-1] + delta;
        
        if (endpoint) v[size-1] = end;
    
    }
    
    return v;
    
}

io::rvec rhistogram(io::rvec const& bins, io::cvec const& data){
    
    io::rvec h(bins.size(),0);
    
    auto const min = bins[0];
    auto const max = bins.back();
    auto const delta=(max-min)/bins.size();
    
    for (std::size_t i=0; i<data.size(); i++){
        auto const index = static_cast<int>((data[i].real()-min)/delta);
        if (index >= bins.size())
            h.back()++;
        else
            h[index]++;
    }
    
    return h;
}

io::rvec rhistogram(io::rvec const& bins, io::rvec const& data){
    
    io::rvec h(bins.size(),0);
    
    auto const min = bins[0];
    auto const max = bins.back();
    auto const delta=(max-min)/bins.size();
    
    for (std::size_t i=0; i<data.size(); i++){
        auto const index = static_cast<int>((data[i]-min)/delta);
        if (index >= bins.size())
            h.back()++;
        else
            h[index]++;
    }
    
    return h;
}

io::rvec ihistogram(io::rvec const& bins, io::cvec const& data){
    
    io::rvec h(bins.size(),0);
    
    auto const min = bins[0];
    auto const max = bins.back();
    auto const delta=(max-min)/bins.size();
    
    for (std::size_t i=0; i<data.size(); i++){
        auto const index = static_cast<int>((data[i].imag()-min)/delta);
        if (index >= bins.size())
            h.back()++;
        else
            h[index]++;
    }
    
    return h;
}

io::rvec ihistogram(io::rvec const& bins, io::rvec const& data){
    return io::rvec(bins.size());
}


template<typename Value>
jsx::value work(jsx::value jParams,double const min, double const max, double const de){
    
    jsx::value jOut;
    
    params::complete_impurity<Value>(jParams);
    
    auto const N = jParams("hybridisation")("matrix").array().size();
    auto const& twoBody = jsx::at<io::Vector<Value>>(jParams("hloc")("two body"));
    auto const& as_twobody = imp::Tensor<Value>(jParams("hloc")("two body"), N);
    
    io::Vector<Value> abba(N*N);
    io::Vector<Value> abab(N*N);
    io::Vector<Value> as_full(N*N*N*N);
    auto it = as_full.begin();
    //Value fingerprint;
    
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++){
            
            int i_abba = i*N*N*N + j*N*N + j*N + i;
            int i_abab = i*N*N*N + j*N*N + i*N + j;
            
            abba[i*N + j] = twoBody[i_abba];
            abab[i*N + j] = twoBody[i_abab];
            

            for (int k=0; k<N; k++)
                for (int l=0; l<N; l++){
                    
                    *it++ = 2.*as_twobody(i,j,k,l);
                    
                }
            
    }
    
    if (min > max)
        throw std::runtime_error("min > max");
    io::rvec bins = linspace( (max-min)/de + 1, min, max);
    
    //jOut["abba"] = abba;
    //jOut["abab"] = abab;
    
    jOut["abba hist"] = rhistogram(bins,abba);
    jOut["abab hist"] = rhistogram(bins,abab);
    
    jOut["full hist"] = rhistogram(bins,twoBody);
    jOut["anti symmetrized full hist"] = rhistogram(bins,as_full);
    
    if (std::is_same<Value,ut::complex>::value)
        jOut["imag full hist"] = ihistogram(bins,twoBody);
    
    jOut["bins"] = bins;
    
    return jOut;
}

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif
    try {
        if(argc != 5) throw std::runtime_error("ctqmc: Wrong number of input parameters: TWOBODY (case name) (min) (max) (de) !");
        
        std::time_t time;  mpi::cout = mpi::cout_mode::one;
        
        mpi::cout << "Start task at " << std::asctime(std::localtime(&(time = std::time(nullptr)))) << std::endl << std::endl;
        
        jsx::value jParams = mpi::read(std::string(argv[1]) + ".json");  params::initialize(jParams); params::complete_worms(jParams);
        
        auto const min = std::stod(argv[2]);
        auto const max = std::stod(argv[3]);
        auto const de = std::stod(argv[4]);
        
        jsx::value jOut;
        if(jParams("complex").boolean()) {
            
            jOut = work<ut::complex>(jParams,min,max,de);
            
        } else {
            
            jOut = work<double>(jParams,min,max,de);
            
        }
        
        mpi::write(jOut,std::string(argv[1]) +".2body.json");
        
    }
    catch (std::exception& exc) {
        std::cerr << exc.what() << " ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        return -1;
    }
    catch (...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
#ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, -2);
#endif
        return -2;
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}












