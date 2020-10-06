driver_dir=./drivers/
evalsim_dir=./evalsim/
ctqmc_dir=./ctqmc/host/
ctqmc_gpu_dir=./ctqmc/device/planar_complex/

include Makefile.in


all:
	+$(MAKE) -C $(evalsim_dir)
	+$(MAKE) -C $(ctqmc_dir)

cpu:
	+$(MAKE) -C $(evalsim_dir)
	+$(MAKE) -C $(ctqmc_dir)

gpu:
	+$(MAKE) -C $(evalsim_dir)
	+$(MAKE) -C $(ctqmc_gpu_dir)

cpu_lib:
	+$(MAKE) -C $(driver_dir) cpu

gpu_lib:
	+$(MAKE) -C $(driver_dir) gpu


clean:
	+$(MAKE) -C $(evalsim_dir) clean
	+$(MAKE) -C $(ctqmc_dir) clean
	+$(MAKE) -C $(ctqmc_gpu_dir) clean

twobody: twobody.C
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -o $@  twobody.C $(LDFLAGS) $(LIBS)
	mv twobody bin/TWOBODY
