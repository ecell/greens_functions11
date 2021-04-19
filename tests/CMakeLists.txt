set(EXE_NAMES
    GreensFunction2DAbsSym
    GreensFunction2DRadAbs
    GreensFunction3DAbsSym
    GreensFunction3DRadAbs
    GreensFunction3DRadInf
    )

foreach(EXE_NAME ${EXE_NAMES})
    add_executable(${EXE_NAME} ${EXE_NAME}.cpp)
    set_target_properties(${EXE_NAME} PROPERTIES
        COMPILE_FLAGS "-O2 -march=native -mtune=native")

    target_link_libraries(${EXE_NAME} greens_functions GSL::gsl GSL::gslcblas)
    if(GF11_BUILD_LIBRARY)
        target_link_libraries(${EXE_NAME} gf11)
    endif()
endforeach(EXE_NAME)
