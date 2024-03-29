# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.14.0)
project (GRIEF_proj C CXX)

# enable_language(CUDA)

include(CheckLanguage)
set(CMAKE_CXX_STANDARD 17)

find_package(PythonLibs REQUIRED)
find_package(OpenCV REQUIRED)
   if(NOT OpenCV_FOUND)
      message(FATAL_ERROR "OpenCV not found.")
   endif()
find_package (Python3 COMPONENTS Interpreter NumPy)
find_package( Eigen3 REQUIRED )

include_directories(
${PROJECT_SOURCE_DIR}
${PROJECT_SOURCE_DIR}/DE
${PROJECT_SOURCE_DIR}/GRIEF_CUDA
${PYTHON_INCLUDE_DIRS}
${EIGEN3_INCLUDE_DIR}
${CUDA_INCLUDE_DIRS}
)

check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
enable_language(CUDA)

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
message("CUDA found!")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
add_library(${PROJECT_NAME} SHARED
DE/DE.cc
GRIEF_CUDA/grief.cu
DE/DE.h
measurements/measurements.h
measurements/measurements.cc
GRIEF_CUDA/grief.h
QS/quicksort.cc
QS/quicksort.h
)


target_link_libraries (${PROJECT_NAME}
${OpenCV_LIBS}
${PYTHON_LIBRARIES}
Python3::NumPy
Eigen3::Eigen

)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (teste ../main.cu)
target_link_libraries(teste ${PROJECT_NAME})
set(CUDA_SEPARABLE_COMPILATION TRUE)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_86,code=sm_86;)
set_target_properties(GRIEF_proj PROPERTIES CUDA_ARCHITECTURES "86")
set_target_properties(teste PROPERTIES CUDA_ARCHITECTURES "86")
set_target_properties( teste PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

else()
message("CUDA not found!")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
add_library(${PROJECT_NAME} SHARED
DE/DE.cc
GRIEF/grief.cc
DE/DE.h
measurements/measurements.h
measurements/measurements.cc
GRIEF/grief.h
QS/quicksort.cc
QS/quicksort.h
)

target_link_libraries (${PROJECT_NAME}
${OpenCV_LIBS}
${PYTHON_LIBRARIES}
Python3::NumPy
Eigen3::Eigen)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (teste ../main.cc)
target_link_libraries(teste ${PROJECT_NAME})

endif()




