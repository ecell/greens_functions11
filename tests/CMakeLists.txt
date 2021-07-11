set(EXE_NAMES
    GreensFunction2DAbsSym_drawTime
    GreensFunction2DAbsSym_drawR

    GreensFunction2DRadAbs_drawTime
    GreensFunction2DRadAbs_drawR

    GreensFunction3DAbsSym_drawTime
    GreensFunction3DAbsSym_drawR

    GreensFunction3DRadAbs_drawTime
    GreensFunction3DRadAbs_drawR
    GreensFunction3DRadAbs_drawTheta

    GreensFunction3DRadInf_drawTime
    GreensFunction3DRadInf_drawR
    GreensFunction3DRadInf_drawTheta
    )

include_directories(${CMAKE_BINARY_DIR})

foreach(EXE_NAME ${EXE_NAMES})
    add_executable(${EXE_NAME} ${EXE_NAME}.cpp)
    set_target_properties(${EXE_NAME} PROPERTIES
        COMPILE_FLAGS "-O2 -march=native -mtune=native")

    target_link_libraries(${EXE_NAME} greens_functions GSL::gsl GSL::gslcblas)
    if(GF11_BUILD_LIBRARY)
        target_link_libraries(${EXE_NAME} gf11)
    endif()
endforeach(EXE_NAME)
