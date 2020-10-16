include ../Makefile.in

include Makefile.objs


BASE_DIR = ../include/
include $(BASE_DIR)Makefile.objs
REL_BASE_OBJS = $(addprefix $(BASE_DIR), $(BASE_OBJS))

MAIN_DIR = ../ctqmc/include/
include $(MAIN_DIR)Makefile.objs
REL_MAIN_OBJS = $(addprefix $(MAIN_DIR), $(MAIN_OBJS))


all:  EVALSIM 
cpu:  EVALSIM
gpu:  EVALSIM

obj: $(EVALSIM_OBJS)

EVALSIM: $(EVALSIM_OBJS)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) $(EVALSIM_OBJS) $(REL_BASE_OBJS) $(REL_MAIN_OBJS) -o $@  Main.C $(LDFLAGS) $(LIBS)
	mv EVALSIM ../bin/.

$(EVALSIM_OBJS):
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -c $(@:.o=.C) -o $@

clean:
	rm -f *.o ../bin/EVALSIM EVALSIM
	find . -type f -name '*.o' -exec rm {} +





 
