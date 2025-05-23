cmake_minimum_required(VERSION 3.1)
FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h" "*.cuh")

if(ENABLE_OMP)
       find_package(OpenMP REQUIRED)
       set(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
       set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_C_FLAGS}")
endif(ENABLE_OMP)
#add_definitions("-Wall -g")

find_package(CUDA QUIET REQUIRED)
find_package(MPI REQUIRED)
# find_package(OpenCV REQUIRED)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR} ${MPI_INCLUDE_PATH})

include_directories(/usr/local/cuda-11.1/targets/x86_64-linux/include)

# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11) # CUDA 使用 C++11 标准
# Pass options to NVCC
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-G -g -O3)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3)
# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -O3 -ffast-math -Wall")
# list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
# list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
# list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_52,code=sm_52)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_52,code=compute_52)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_62,code=sm_62)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_80,code=sm_80)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_86,code=sm_86)
#-G -g -O0 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -c  -o

cuda_add_executable(Markerfree ${SOURCES})
target_link_libraries(Markerfree cufft ${OpenCV_LIBS})

# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")

target_link_libraries(Markerfree ${MPI_LIBRARIES} method_library mrc)