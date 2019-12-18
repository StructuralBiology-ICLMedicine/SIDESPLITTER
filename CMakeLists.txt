cmake_minimum_required(VERSION 3.1)

project(SIDESPLITTER C)

# Set default build type to Release (implies -O3)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Find pthread library
find_package(Threads REQUIRED)

# Find FFTW3
find_package(FFTW3 CONFIG REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})
find_library(FFTW3 fftw3)
find_library(FFTW3_THREADS fftw3_threads)

# Find all source files
file(GLOB SOURCES "*.c")

# Compile and link the executable
add_executable(sidesplitter ${SOURCES})
set_property(TARGET sidesplitter PROPERTY C_STANDARD 99)
target_link_libraries(sidesplitter m ${FFTW3} ${FFTW3_THREADS} ${CMAKE_THREAD_LIBS_INIT})

# Install the executable into the bin directory
# (This is done relative to the CMAKE_INSTALL_PREFIX and is a sensible default)
install(TARGETS sidesplitter DESTINATION bin)

# Install the wrapper script
install(PROGRAMS sidesplitter_wrapper.sh DESTINATION bin)

