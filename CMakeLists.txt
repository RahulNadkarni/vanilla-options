cmake_minimum_required(VERSION 3.10)
project(VanillaOptions)

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(vanilla-options
    src/main.cpp
    src/Utils.cpp
    src/BlackScholes.cpp
    src/Options.cpp
    src/MonteCarlo.cpp
    src/EulerMonteCarlo.cpp 
)

add_executable(consistency-tests 
    tests/ConsistencyTests.cpp
    src/Utils.cpp 
    src/BlackScholes.cpp)