include ../../Makefile.in

all: CTQMC
cpu: CTQMC

CTQMC:  ctqmc.C $(HEADERS_IS)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -o $@ ctqmc.C $(LDFLAGS) $(LIBS)
	mv CTQMC ../../bin/.

clean:
	rm -f *.o ../../bin/CTQMC CTQMC
	









 
