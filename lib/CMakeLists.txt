set (S2D_OBJS utils.cpp hmc.cpp inverters.cpp dirac_op.cpp measurements.cpp blas.cpp)

add_library(schwinger2d STATIC ${S2D_OBJS})
target_include_directories(schwinger2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)

#add_library(quda STATIC $<TARGET_OBJECTS:quda_cpp> ${QUDA_CU_OBJS})
