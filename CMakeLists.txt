cmake_minimum_required(VERSION 3.10)
project(VanillaOptions)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

set(SOURCES
    src/BlackScholes.cpp
    src/Options.cpp 
    src/MonteCarlo.cpp 
    src/EulerMonteCarlo.cpp 
    src/Utils.cpp
    )

add_executable(consistency_tests
    tests/ConsistencyTests.cpp
    ${SOURCES}
    )

add_executable(montecarlo_tests 
    tests/MonteCarloComp.cpp
    ${SOURCES}
    )
    