include ../include/msource/MakefileHomeMPI.in

HEADERS =  

all:  EVALSIM 

EVALSIM:  Evalsim.C $(HEADERS)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -o $@  Evalsim.C $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o EVALSIM 





 
