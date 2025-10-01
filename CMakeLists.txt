cmake_minimum_required(VERSION 3.10)
project(VanillaOptions)

set(CMAKE_CXX_STANDARD 17)

# Enable optimization flags for maximum performance
# -O3: Maximum optimization level
# -march=native: Optimize for the specific CPU architecture
# -DNDEBUG: Disable debug assertions in release mode
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_BUILD_TYPE Release)

# Find and enable OpenMP for parallel Monte Carlo simulations
# OpenMP provides multi-threading capabilities for the optimized Monte Carlo methods
if(APPLE)
    # For macOS, manually set OpenMP paths (libomp from Homebrew)
    # This enables parallel processing in Monte Carlo simulations
    set(OpenMP_C_FLAGS "-Xclang -fopenmp -I/usr/local/opt/libomp/include")
    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I/usr/local/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
    set(OpenMP_CXX_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
    set(OpenMP_CXX_FOUND TRUE)
    set(OpenMP_C_FOUND TRUE)
else()
    find_package(OpenMP REQUIRED)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

# Main executable
add_executable(vanilla_options
    src/main.cpp
    src/BlackScholes.cpp
    src/Options.cpp
    src/MonteCarlo.cpp
    src/EulerMonteCarlo.cpp
    src/Utils.cpp
    src/Greeks.cpp
)

# Link OpenMP to the main executable for parallel Monte Carlo simulations
# This enables the parallel processing optimization methods
if(APPLE)
    target_link_libraries(vanilla_options "/usr/local/opt/libomp/lib/libomp.dylib")
    target_compile_options(vanilla_options PRIVATE "-Xclang" "-fopenmp" "-I/usr/local/opt/libomp/include")
else()
    target_link_libraries(vanilla_options OpenMP::OpenMP_CXX)
endif()

# Consistency tests
add_executable(consistency_tests 
    tests/ConsistencyTests.cpp
    src/Utils.cpp 
    src/BlackScholes.cpp
    src/Options.cpp
    src/Greeks.cpp
)

# Link OpenMP to consistency tests for parallel Monte Carlo validation
if(APPLE)
    target_link_libraries(consistency_tests "/usr/local/opt/libomp/lib/libomp.dylib")
    target_compile_options(consistency_tests PRIVATE "-Xclang" "-fopenmp" "-I/usr/local/opt/libomp/include")
else()
    target_link_libraries(consistency_tests OpenMP::OpenMP_CXX)
endif()
