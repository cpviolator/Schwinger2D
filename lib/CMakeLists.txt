set (S2D_OBJS utils.cpp hmc.cpp inverters.cpp dirac_op.cpp measurements.cpp blas.cpp iram.cpp)

add_library(schwinger2d STATIC ${S2D_OBJS})
target_include_directories(schwinger2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)
include_directories(${EIGEN_INCLUDE_DIRS})
add_dependencies(schwinger2d Eigen)
add_dependencies(schwinger2d HDF5)


if(ENABLE_OPENMP)
  target_link_libraries(schwinger2d PUBLIC OpenMP::OpenMP_CXX)
endif()
