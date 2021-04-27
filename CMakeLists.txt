cmake_minimum_required(VERSION 3.19)

project(ComCTQMC VERSION 1.0.0)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

include(CheckLanguage)


##################################################################


enable_language(CXX)

check_language(CUDA)
IF(CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
	add_definitions(-DUSE_CUDA_SINGLE)

	include_directories(${CUDA_INCLUDE_PATH})
	add_link_options(${CUDA_LINKER_FLAGS})
ELSE()
	message(STATUS "No CUDA support detected")
ENDIF()

find_package(MPI)
IF(MPI_FOUND)
	add_definitions(-DHAVE_MPI)
	include_directories(${MPI_INCLUDE_PATH})
	add_link_options(${MPI_LINKER_FLAGS})
ENDIF()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_link_options(${BLAS_LINKER_FLAGS})

IF(NOT CMAKE_BUILD_TYPE EQUAL "DEBUG")
	add_definitions(-DNDEBUG)
ENDIF()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_definitions(-DCMAKE)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

##################################################################

add_library(ctqmc SHARED)
target_sources(ctqmc PUBLIC lib/libCTQMC.C)
target_link_libraries(ctqmc
    ${MPI_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

add_executable(ctqmc_x)
if (CMAKE_CUDA_COMPILER)
	target_sources(ctqmc_x PUBLIC ctqmc/device/planar/complex/Main.C)
ELSE()
	target_sources(ctqmc_x PUBLIC ctqmc/host/Main.C)
ENDIF()
target_link_libraries(ctqmc_x ctqmc
    ${MPI_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

add_executable(evalsim_x)
target_sources(evalsim_x PUBLIC evalsim/Main.C)
target_link_libraries(evalsim_x ctqmc
    ${MPI_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

##################################################################

add_subdirectory(include)
add_subdirectory(evalsim)
add_subdirectory(ctqmc)


