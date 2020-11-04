include ../../Makefile.in
include Makefile.objs

all : $(MAIN_OBJS)

$(MAIN_OBJS):
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -c $(@:.o=.C) -o $@

clean:
	find . -type f -name '*.o' -exec rm {} +

	









 
