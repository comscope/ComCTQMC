#ifndef INCLUDE_FOURIERTRANSFORM
#define INCLUDE_FOURIERTRANSFORM

/* Translated From Portobello's Maxent */

namespace ft{

    template <typename T>
    T transformF(T const t, std::vector<std::complex<T>> const& Gm, std::vector<T> const& om, T const ah, T const beta){
        T dsum=0.;
        for (int im=0; im < Gm.size(); im++){
            dsum += std::cos(om[im]*t)*std::real(Gm[im]) + std::sin(om[im]*t) * (std::imag(Gm[im])+ah/om[im]);
        }
        return 2*dsum/beta-0.5*ah;
    }

    template <typename T>
    T transformB(T const t, std::vector<std::complex<T>> const& Gm, std::vector<T> const& om, T const beta){
        double dsum=0.;
        for (int im=1; im < Gm.size(); im++)
            dsum += std::cos(om[im]*t)*std::real(Gm[im]);
        dsum = dsum + 0.5*std::real(Gm[1]);
        return 2*dsum/beta;
    }

    template <typename T>
    T findHighFrequency(std::vector<std::complex<T>> const& Gm, std::vector<T> const& om, int const Nf){
        T S=0.;T Sx=0.;T Sy=0.;T Sxx=0.;T Sxy=0;
        for (int j = om.size() - Nf; j < om.size(); j++){
            auto x = om[j];
            auto y = Gm[j].imag() * x;
            auto x2 = x*x;
            
            Sy += y;
            Sx += 1./x2;
            Sxx += 1./(x2*x2);
            Sxy += y/x2;
            S += 1;
        }
            
        auto dd = S*Sxx-Sx*Sx;
        auto a = (Sxx*Sy-Sx*Sxy)/dd;
        auto ah = -a;
            
        if (std::abs(ah-1.0)<1e-3)
            ah=1.0;
        
        return ah;
    }

    template <typename T>
    std::vector<std::complex<T>> inverseFourier(std::vector<std::complex<T>> const& Gm, std::vector<T> const& om, std::vector<T> tau, T const beta, int const Nf, std::string const& stat = "fermi"){
        
        std::vector<std::complex<T>> Gtau(tau.size());
       
        auto df = Gm.back().real() * om.back()/M_PI;
        
        if (stat=="fermi"){
            auto ah = findHighFrequency(Gm,om,Nf);
            
            for (int it=0; it < tau.size(); it++){
                auto const t = tau[it];
                Gtau[it] = transformF(t,Gm,om,ah,beta);
            }
            Gtau[0] += df;
            Gtau.back() -= df;
        } else {

            for (int it=0; it < tau.size(); it++){
                auto const t = tau[it];
                Gtau[it] = transformB(t,Gm,om,beta);
            }
            Gtau[0] += df;
            Gtau.back() += df;
        }
        
        return Gtau;

    }
    
}


#endif
