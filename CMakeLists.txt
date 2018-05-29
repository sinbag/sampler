cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 OLD)
#list(APPEND CMAKE_MODULE_PATH "cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/")
##### BUILD ############################################################
project(sampler)

option(CUDA_ENABLED "If machine supports CUDA / NVIDIA Graphics Card then run
cmake -DCUDA_ENABLED=ON .." OFF)

if(CUDA_ENABLED)
add_definitions(-DUSE_CUDA_ENABLED)
endif(CUDA_ENABLED)

if(CUDA_ENABLED)
find_package(CUDA REQUIRED)
endif(CUDA_ENABLED)


#set(EXECUTABLE_OUTPUT_PATH ./)

#find_package(Boost 1.42.0 COMPONENTS program_options REQUIRED)
find_package(Boost 1.42.0 COMPONENTS system filesystem date_time program_options REQUIRED)
#find_package(Qt5 COMPONENTS QtCore QtGui QtOpenGl REQUIRED)
#find_package(OpenGL REQUIRED)
find_package(OpenEXR)
find_package(IlmBase)
find_package(TBB)

INCLUDE(FindPkgConfig)

PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
PKG_SEARCH_MODULE(SDL2IMAGE REQUIRED SDL2_image>=2.0.0)

#include(${QT_USE_FILE})

include_directories(
        ${PROJECT_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR}
#        ${QT_QTOPENGL_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR}
        ${PROJECT_BINARY_DIR}
        ${SDL2_INCLUDE_DIR}
        ${SDL2IMAGE_INCLUDE_DIR}
  ${OPENEXR_INCLUDE_DIR}
  ${ILMBASE_INCLUDE_DIR}


)

link_directories(
        ${Boost_LIBRARY_DIRS}
        ${OPENGL_LIBRARY_DIRS}
)



######## External Data for production/barycenter/lut files######
##
## Data files will be located in the {build-folder}/data/lut/
##
include(ExternalData)
set (ExternalData_TIMEOUT_ABSOLUTE 0)
list(APPEND ExternalData_URL_TEMPLATES
  "http://liris.cnrs.fr/polyhex/data/%(algo)-%(hash)")
ExternalData_Expand_Arguments(
  GetData
  OUT-data
  DATA{data/lut/production_rules.dat}
  DATA{data/lut/barycenters.dat}
  )
MESSAGE(STATUS "PATH to the data files: ${OUT-data}")
ExternalData_Add_Target(GetData)
ExternalData_Expand_Arguments(
  GetLUT
  OUT-LUT
  DATA{data/lut/offsets_bnot.dat}
  )
MESSAGE(STATUS "PATH to the LUT file: ${OUT-LUT}")
ExternalData_Add_Target(GetLUT)
set_target_properties(GetLUT PROPERTIES EXCLUDE_FROM_ALL "TRUE")
##########"

add_definitions(-std=gnu++0x)

if(CUDA_FOUND)

include_directories(
#CUDA files
${CMAKE_CURRENT_SOURCE_DIR}/src/cuda
)

file ( GLOB SDNN_CUDA_FILES
src/cuda/*.cu
src/cuda/*.h
src/cuda/*.cpp
)

# compared to class settings, we let NVidia's FindCUDA CMake detect
# whether to build x64.  We tell it to support most devices.

set(CUDA_NVCC_FLAGS "-gencode;arch=compute_61,code=sm_61;")

# add -Wextra compiler flag for gcc compilations
if (UNIX)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-Xcompiler -Wextra")
#set(CMAKE_CXX_FLAGS "-stdlib=libstdc++")
endif (UNIX)

# add debugging to CUDA NVCC flags.  For NVidia's NSight tools.
set(CUDA_NVCC_FLAGS_DEBUG ${CUDA_NVCC_FLAGS_DEBUG} "-G")


add_executable(fourier-analyzer-cuda ${PROJECT_SOURCE_DIR}/src/main/cft-cuda-demo.cxx ${other-src} ${sampler-core-src} ${SDNN_CUDA_FILES})
target_link_libraries(
fourier-analyzer-cuda
${Boost_LIBRARIES}
${SDL2_LIBRARIES}
${SDL2IMAGE_LIBRARIES}
${OPENEXR_LIBRARIES}
${ILMBASE_LIBRARIES}
${TBB_LIBRARIES}
)

else(CUDA_FOUND)
message("CUDA is not installed on this system.")
endif(CUDA_FOUND)


file(GLOB sampler-core-src
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB other-src
  ${PROJECT_SOURCE_DIR}/src/io/*.cpp
  ${PROJECT_SOURCE_DIR}/src/core/*.cpp
  ${PROJECT_SOURCE_DIR}/src/fourier-analysis/*.cpp
  ${PROJECT_SOURCE_DIR}/src/util-samples/*.cpp
  ${PROJECT_SOURCE_DIR}/src/statistics/*.cpp
  ${PROJECT_SOURCE_DIR}/src/test-functions/*.cpp
)

file(GLOB sampler-gui-src
${PROJECT_SOURCE_DIR}/src/gui/*.cpp
)
file(GLOB sampler-gui-hdr
${PROJECT_SOURCE_DIR}/src/gui/*.hpp
)

add_executable(sampler-cli-demo1 ${PROJECT_SOURCE_DIR}/src/main/cli-demo1.cxx ${other-src} ${sampler-core-src})
target_link_libraries(
        sampler-cli-demo1
        ${Boost_LIBRARIES}
        ${SDL2_LIBRARIES}
        ${SDL2IMAGE_LIBRARIES}
${OPENEXR_LIBRARIES}
${ILMBASE_LIBRARIES}
${TBB_LIBRARIES}
)

add_executable(sampler-cli ${PROJECT_SOURCE_DIR}/src/main/cli.cpp ${sampler-core-src})
target_link_libraries(
        sampler-cli
        ${Boost_LIBRARIES}
        ${SDL2_LIBRARIES}
        ${SDL2IMAGE_LIBRARIES}
)
add_executable(fourier-analysis ${PROJECT_SOURCE_DIR}/src/main/cli-fourier-analysis-demo1.cxx ${other-src} ${sampler-core-src})
target_link_libraries(
        fourier-analysis
        ${Boost_LIBRARIES}
        ${SDL2_LIBRARIES}
        ${SDL2IMAGE_LIBRARIES}
        ${OPENEXR_LIBRARIES}
        ${ILMBASE_LIBRARIES}
        ${TBB_LIBRARIES}
)

add_executable(variance-analysis ${PROJECT_SOURCE_DIR}/src/main/variance-analysis-demo1.cxx ${other-src} ${sampler-core-src})
target_link_libraries(
        variance-analysis
        ${Boost_LIBRARIES}
        ${SDL2_LIBRARIES}
        ${SDL2IMAGE_LIBRARIES}
        ${OPENEXR_LIBRARIES}
        ${ILMBASE_LIBRARIES}
        ${TBB_LIBRARIES}
)

#qt4_wrap_cpp(sampler-gui-moc ${sampler-gui-hdr})
#add_executable(sampler-gui src/main/gui.cpp ${sampler-gui-moc} ${sampler-gui-src} ${sampler-core-src})
#target_link_libraries(
#        sampler-gui
#        ${Boost_LIBRARIES}
#        ${QT_LIBRARIES}
#        ${OPENGL_LIBRARY}
#        ${SDL2_LIBRARIES}
#        ${SDL2IMAGE_LIBRARIES}
#)

add_executable(sampler-test src/main/test.cpp ${sampler-core-src})
target_link_libraries(
        sampler-test
        ${Boost_LIBRARIES}
        ${SDL2_LIBRARIES}
        ${SDL2IMAGE_LIBRARIES}
)

install(
        PROGRAMS
#        ./bin/sampler-cli
#        ./bin/sampler-gui
        DESTINATION bin
)

#------------------------------------------------------------------------------
# Configuration of the Config.h
#------------------------------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/src/config.h.in
  ${PROJECT_BINARY_DIR}/config.h)
