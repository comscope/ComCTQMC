include ../Makefile.in
include Makefile.objs

all : $(BASE_OBJS)

$(BASE_OBJS):
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -c $(@:.o=.C) -o $@

clean:	
	find . -type f -name '*.o' -exec rm {} +
	









 
