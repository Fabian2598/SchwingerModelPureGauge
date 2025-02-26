cmake_minimum_required(VERSION 3.29.5)
set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_COMPILER "C:/msys64/ucrt64/bin/gcc.exe")
set(CMAKE_CXX_COMPILER "C:/msys64/ucrt64/bin/g++.exe")

project(Schwinger_model VERSION 1.0)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/statistics.cpp
    src/gauge_conf.cpp
    src/variables.cpp
    src/metropolis.cpp
    src/main.cpp   
    src/hmc.cpp
)

set(NS "100") 
set(NT "100") 
configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)

# Add an executable with the above sources
add_executable("SM_${NS}x${NT}" ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_options("SM_${NS}x${NT}" PRIVATE
			-O3
)

#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../