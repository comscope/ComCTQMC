add_subdirectory(include)
add_subdirectory(host)
if (CMAKE_CUDA_COMPILER)
	add_subdirectory(device/planar_complex)
ENDIF()

