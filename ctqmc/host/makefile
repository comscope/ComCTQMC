include ../../Makefile.in

all: CTQMC
cpu: CTQMC

SUB_OBJS_FILES = \
mpi.o \
evalsim.o \
shared.o 

OBJS_FILES = \
dhMonteCarlo.o \
zhMontecarlo.o

OBJS_DIR = obj/

OBJS = $(addprefix $(OBJS_DIR), $(OBJS_FILES))
SUB_OBJS = $(addprefix $(OBJS_DIR), $(SUB_OBJS_FILES))

CTQMC : $(SUB_OBJS_FILES) $(OBJS_FILES)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS)  $(SUB_OBJS) $(OBJS) -o $@ ctqmc.C  $(LDFLAGS) $(LIBS)
	mv CTQMC ../../bin/.

$(SUB_OBJS_FILES):
	@mkdir -p $(OBJS_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -c $(@:.o=.C) -o $(OBJS_DIR)$@

$(OBJS_FILES):
	@mkdir -p $(OBJS_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -c $(@:.o=.C) -o $(OBJS_DIR)$@

clean:
	rm -f *.o ../../bin/CTQMC CTQMC
	









 
