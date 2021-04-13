set(EXE_NAMES
    GreensFunction2DAbsSym
    GreensFunction2DRadAbs
    GreensFunction3DAbsSym
    GreensFunction3DRadAbs
    GreensFunction3DRadInf
    )

find_package(Boost 1.67.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

foreach(EXE_NAME ${EXE_NAMES})
    add_executable(${EXE_NAME} ${EXE_NAME}.cpp)
    set_target_properties(${EXE_NAME} PROPERTIES
        COMPILE_FLAGS "-O2 -march=native -mtune=native")

    target_link_libraries(${EXE_NAME} greens_functions GSL::gsl GSL::gslcblas)
endforeach(EXE_NAME)
