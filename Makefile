evalsim_dir=./evalsim/
ctqmc_dir=./ctqmc/host/
ctqmc_gpu_dir=./ctqmc/device/planar_complex/

clean:
	+$(MAKE) -C $(evalsim_dir) clean
	+$(MAKE) -C $(ctqmc_dir) clean
	+$(MAKE) -C $(ctqmc_gpu_dir) clean

all:
	+$(MAKE) -C $(evalsim_dir)
	+$(MAKE) -C $(ctqmc_gpu_dir)

cpu:
	+$(MAKE) -C $(evalsim_dir)
	+$(MAKE) -C $(ctqmc_dir)

gpu:
	+$(MAKE) -C $(ctqmc_dir)
	+$(MAKE) -C $(ctqmc_gpu_dir)

