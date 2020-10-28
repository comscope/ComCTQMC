#if we haven't already set the necessary options, try loading them from Makefile.in
ifndef CXX_MPI
include Makefile.in
endif

ifndef EXE_DIR
EXE_DIR = ./bin/
endif

ifndef LIB_DIR
LIB_DIR = ./lib/
endif 


EVALSIM_DIR=./evalsim/
HOST_DIR=./ctqmc/host/
GPU_DIR=./ctqmc/device/planar_complex/
BASE_DIR=./include/
MAIN_DIR=./ctqmc/include/

include $(EVALSIM_DIR)makefile.objs
include $(HOST_DIR)makefile.objs
include $(GPU_DIR)makefile.objs
include $(BASE_DIR)makefile.objs
include $(MAIN_DIR)makefile.objs

MAIN_OBJS := $(addprefix $(MAIN_DIR), $(MAIN_OBJS))
BASE_OBJS := $(addprefix $(BASE_DIR), $(BASE_OBJS))
EVALSIM_OBJS := $(addprefix $(EVALSIM_DIR), $(EVALSIM_OBJS))

HOST_OBJS := $(addprefix $(HOST_DIR), $(HOST_OBJS))

GPU_OBJS := $(addprefix $(GPU_DIR), $(GPU_OBJS))
GPU_CUDA_OBJS := $(addprefix $(GPU_DIR), $(GPU_CUDA_OBJS))

all : cpu evalsim

cpu : CTQMC_x

gpu : CTQMC_gx

evalsim : EVALSIM_x

lib : CXXFLAGS += -fPIC
lib : libCTQMC.so

%.o : %.c
	$(CXX_MPI) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

libCTQMC.so : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS)
	@mkdir -p $(LIB_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -shared $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) -o $(LIB_DIR)/$@ ./lib/libCTQMC.C $(LFLAGS) $(LIBS)

EVALSIM_x : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS)
	@mkdir -p $(EXE_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) -o $(EXE_DIR)$@ $(EVALSIM_DIR)Main.C $(LDFLAGS) $(LIBS)
	cp $(EXE_DIR)$@ $(EXE_DIR)EVALSIM

CTQMC_x : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) 
	@mkdir -p $(EXE_DIR)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) -o $(EXE_DIR)$@ $(HOST_DIR)Main.C $(LDFLAGS) $(LIBS)
	cp $(EXE_DIR)$@ $(EXE_DIR)CTQMC

CTQMC_gx : CXXFLAGS += -DMAKE_GPU_ENABLED
CTQMC_gx : $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) $(GPU_OBJS) $(GPU_CUDA_OBJS)
	@mkdir -p $(EXE_DIR)
	$(CXX_MPI) -o $(EXE_DIR)$@  $(BASE_OBJS) $(MAIN_OBJS) $(EVALSIM_OBJS) $(HOST_OBJS) $(GPU_OBJS) $(GPU_CUDA_OBJS) $(GPU_CUDA_OBJS)bj.o $(GPU_DIR)Main.C -lcudart -lcudadevrt $(LDFLAGS) $(LIBS)
	cp $(EXE_DIR)$@ $(EXE_DIR)CTQMC

$(GPU_CUDA_OBJS) :
	$(NVCC) $(NVCCFLAGS) $(CUDA_CPPFLAGS) $(CUTLASS_CPPFLAGS) -dc $(@:.o=.cu) -o $@bj.o
	$(NVCC) -arch=$(GPU_ARCH) -dlink -lcublas_device $@bj.o -o $@
	
clean : 
	rm -f $(EXE_DIR)/EVALSIM* $(EXE_DIR)/CTQMC* $(LIB_DIR)/*.so
	find . -type f -name '*.o' -exec rm {} +






