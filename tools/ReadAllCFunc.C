#include <algorithm>
#include <cmath>

#include "../include/JsonX.h"
#include "../include/mpi/Utilities.h"
#include "../include/io/Vector.h"

int main(int argc, char** argv)
{
	try {		
		if(argc < 2)
			throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
		

        jsx::value jObservables = mpi::read(std::string(argv[1]) + ".json");

        for(int i = 2; i < argc; ++i) {
            jsx::value temp = jObservables(argv[i]);
            jObservables = std::move(temp);
        }

        std::size_t const N = jsx::at<io::cvec>(jObservables.object().begin()->second).size();

        std::ofstream file((std::string(argv[argc - 1]) + ".dat").c_str());
        
        
        int counter = 2;
        for(auto& jFunction : jObservables.object()) {
            std::cout << jFunction.first << " : " << counter << std::endl;
            counter += 2;
        }
        
		for(std::size_t n = 0; n < N; ++n) {
			file << n;
            
            for(auto& jFunction : jObservables.object())
				file << " " << jsx::at<io::cvec>(jFunction.second)[n].real() << " " << jsx::at<io::cvec>(jFunction.second)[n].imag() << " ";
            
			file << std::endl;
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












