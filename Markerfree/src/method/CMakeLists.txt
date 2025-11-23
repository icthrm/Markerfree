
file(GLOB METHOD_SOURCES "*.cpp" "*.cu" "*.h" "*.cuh")
message(STATUS "METHOD_SOURCES: ${METHOD_SOURCES}")

if(ENABLE_OMP)
       find_package(OpenMP REQUIRED)
       set(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
       set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_C_FLAGS}")
endif(ENABLE_OMP)

find_package(CUDA QUIET REQUIRED)
find_package(MPI REQUIRED)
# find_package(OpenCV REQUIRED)

# set(ITK_DIR /home/liuzh/ITK/build)
# find_package(ITK REQUIRED)
# include(${ITK_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR} ${MPI_INCLUDE_PATH})

include_directories(/usr/local/cuda-11.1/targets/x86_64-linux/include)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11) # CUDA 使用 C++11 标准
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})
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

cuda_add_library(method_library STATIC ${METHOD_SOURCES})
target_link_libraries(method_library cufft)
target_link_libraries(method_library ${MPI_LIBRARIES} mrc)

# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")


