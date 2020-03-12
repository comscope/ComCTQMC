#include <algorithm>
#include "../include/io/Vector.h"
#include "../include/mpi/Utilities.h"

int main(int argc, char** argv)
{
	try {		
		if(argc != 7)
			throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
        
        /*
        double const beta = std::atof(argv[1]);
        int const j1 = atoi(argv[2]);
        int const j2 = atoi(argv[3]);
        
        std::vector<std::vector<double>> mj{{-2.5, -1.5, -0.5, 0.5, 1.5, 2.5}, {-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5}};
        std::vector<std::vector<int>> index{{0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11, 12, 13}};
        
        std::vector<std::string> entries;
        
        jsx::value jObservables = mpi::read("params.obs.json");
        
        std::size_t min = jsx::at<io::rvec>(jObservables("occupation-susceptibility")("0_0")("function")).size();
        
        std::ofstream file("susc.dat");
        
        for(std::size_t n = 0; n < min; ++n) {
            file << 2*n*M_PI/beta;
            
            double sum = .0;
            
            for(std::size_t mj1 = 0; mj1 < mj[j1].size(); ++mj1)
                for(std::size_t mj2 = 0; mj2 < mj[j2].size(); ++mj2)
                    sum += mj[j1][mj1]*mj[j2][mj2]*jsx::at<io::rvec>(jObservables("occupation-susceptibility")(std::to_string(index[j1][mj1]) + "_" + std::to_string(index[j2][mj2]))("function"))[n];
            
            file << " " << sum << std::endl;
        }
        
        file.close();
        */
        
        
        double const beta = std::atof(argv[1]);
        std::string name(argv[2]);
        int const sA = atoi(argv[3]);
        int const eA = atoi(argv[4]);
        int const sB = atoi(argv[5]);
        int const eB = atoi(argv[6]);
        
        std::vector<std::string> entries;

        jsx::value jObservables = mpi::read("params.obs.json");
        
        std::size_t min = jsx::at<io::rvec>(jObservables("occupation-susceptibility-" + name)("0_0")("function")).size();
        
        std::ofstream file("sum_" + name + std::to_string(sA) + "-" + std::to_string(eA) + "_" + std::to_string(sB) + "-" + std::to_string(eB) + ".dat");
        
        for(std::size_t n = 0; n < min; ++n) {
            file << 2*n*M_PI/beta;
            
            double sum = .0;
            
            for(int i = sA; i < eA; ++i)
                for(int j = sB; j < eB; ++j)
                    sum += jsx::at<io::rvec>(jObservables("occupation-susceptibility-" + name)(std::to_string(i) + "_" + std::to_string(j))("function"))[n];

            file << " " << sum << std::endl;
        }
        
        file.close();
	}
	catch(std::exception& exc) {
		std::cerr << exc.what() << "\n";
		return -1;
	}
	catch(...) {
		std::cerr << "Fatal Error: Unknown Exception!\n";
		return -2;
	}
	return 0;
}












