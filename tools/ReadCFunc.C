#include <algorithm>
#include <cmath>

#include "../include/JsonX.h"
#include "../include/mpi/Utilities.h"
#include "../include/io/Vector.h"

int main(int argc, char** argv)
{
	try {		
		if(argc < 4)
			throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
		
		std::string name = argv[1];
		double const beta = std::atof(argv[2]);
		
		std::vector<std::string> entries;
		
		for(int i = 3; i < argc; ++i) entries.push_back(argv[i]);

        jsx::value jObservables = mpi::read("params.obs.json")("partition");

		std::size_t min = jsx::at<io::cvec>(jObservables(name)(entries[0])("function")).size();
		for(std::size_t i = 0; i < entries.size(); ++i) 
			min = std::min(min, jsx::at<io::cvec>(jObservables(name)(entries[i])("function")).size());

		std::ofstream file((name + ".dat").c_str());
		
		for(std::size_t n = 0; n < min; ++n) { 
			file << (2*n + 1)*M_PI/beta;
			for(std::size_t i = 0; i < entries.size(); ++i) 
				file << " "
                     << jsx::at<io::cvec>(jObservables(name)(entries[i])("function"))[n].real()
                     << " "
                     << jsx::at<io::cvec>(jObservables(name)(entries[i])("function"))[n].imag();
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












