cmake_minimum_required(VERSION 3.8)
project(gf11 LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/extlib)
include_directories(${PROJECT_SOURCE_DIR}/extlib/greens_functions/build)

find_package(Boost 1.67.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

option(GF11_BUILD_LIBRARY "build library" OFF)

if(GF11_BUILD_LIBRARY)
    message(STATUS "building library")
    add_library(gf11 STATIC
        include/gf11/config.hpp
        include/gf11/factorial.hpp
        include/gf11/find_root.hpp
        include/gf11/gf_math.hpp
        include/gf11/GreensFunction2DAbsSym.hpp
        include/gf11/GreensFunction2DAbsSym.cpp
        include/gf11/GreensFunction2DRadAbs.hpp
        include/gf11/GreensFunction2DRadAbs.cpp
        include/gf11/GreensFunction3DAbsSym.hpp
        include/gf11/GreensFunction3DAbsSym.cpp
        include/gf11/GreensFunction3DRadAbs.hpp
        include/gf11/GreensFunction3DRadInf.hpp
        include/gf11/index_sequence.hpp
        include/gf11/ipv_event_kind.hpp
        include/gf11/SphericalBesselGenerator.hpp
        include/gf11/sumup.hpp
        include/gf11/throw_exception.hpp
        include/gf11/tolerance.hpp
        include/gf11/type_traits.hpp
        )
else()
    message(STATUS "defining GF11_HEADER_ONLY")
    add_definitions("-DGF11_HEADER_ONLY")
endif()

add_subdirectory(extlib)
add_subdirectory(src)
