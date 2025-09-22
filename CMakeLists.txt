cmake_minimum_required(VERSION 3.10)
project(VanillaOptions)

set(CMAKE_CXX_STANDARD 17)

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

# Consistency tests
add_executable(consistency_tests 
    tests/ConsistencyTests.cpp
    src/Utils.cpp 
    src/BlackScholes.cpp
    src/Options.cpp
    src/Greeks.cpp
)
