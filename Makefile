#if we haven't already set the necessary options, try loading them from Makefile.in
ifndef CXX_MPI
include Makefile.in
endif

EVALSIM_DIR=./evalsim/
HOST_DIR=./ctqmc/host/
GPU_DIR=./ctqmc/device/planar_complex/
BASE_DIR=./include/
MAIN_DIR=./ctqmc/include/

include $(EVALSIM_DIR)Makefile.objs
include $(HOST_DIR)Makefile.objs
include $(GPU_DIR)Makefile.objs
include $(BASE_DIR)Makefile.objs
include $(MAIN_DIR)Makefile.objs

MAIN_OBJS := $(addprefix $(MAIN_DIR), $(MAIN_OBJS))
BASE_OBJS := $(addprefix $(BASE_DIR), $(BASE_OBJS))
EVALSIM_OBJS := $(addprefix $(EVALSIM_DIR), $(EVALSIM_OBJS))

HOST_OBJS := $(addprefix $(HOST_DIR), $(HOST_OBJS))

GPU_C_OBJS := $(addprefix $(GPU_DIR), $(GPU_C_OBJS))
GPU_CUDA_OBJS := $(addprefix $(GPU_DIR), $(GPU_C_OBJS))

all : cpu evalsim

cpu : CTQMC_x

gpu : ADD_FLAG CTQMC_gx

evalsim : EVALSIM_x

lib : libCTQMC.so

ADD_FLAG :
	CXXFLAGS = $(CXXFLAGS) -DMAKE_GPU_ENABLED

%.o : %.c
	$(CXX_MPI) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

libCTQMC.so : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS)
	@mkdir -p $(LIB_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -shared $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) -o $(LIB_DIR)/$@ ./lib/libCTQMC.C $(LFLAGS) $(LIBS)

EVALSIM_x : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS)
	@mkdir -p $(EXE_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) -o $(EXE_DIR)$@ $(EVALSIM_DIR)Main.C $(LDFLAGS) $(LIBS)

CTQMC_x : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) 
	@mkdir -p $(EXE_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) -o $(EXE_DIR)$@ $(HOST_DIR)Main.C $(LDFLAGS) $(LIBS)

CTQMC_gx : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(GPU_C_OBJS) $(GPU_CUDA_OBJS)
	@mkdir -p $(EXE_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(GPU_C_OBJS) $(GPU_CUDA_OBJS) -o $(EXE_DIR)$@ $(GPU_DIR)Main.C $(LDFLAGS) $(LIBS)

$(GPU_CUDA OBJS) :
	$(NVCC) $(NVCCFLAGS_lib) $(CUDA_CPPFLAGS) $(CUTLASS_CPPFLAGS) -dc $(@:.o=.cu) -o $@bj.o
	$(NVCC) $(GPU_ARCH) -dlink $@bj.o -o $@ 
	
clean : 
	rm -f $(EXE_DIR)/EVALSIM_x $(EXE_DIR)/CTQMC_x $(EXE_DIR)/CTQMC_gx $(LIB-DIR)/*.so
	find . -type f -name '*.o' -exec rm {} +






