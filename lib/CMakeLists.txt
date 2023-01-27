set (SCHWINGER_2D_OBJS utils.cpp hmc.cpp inverters.cpp dirac_op.cpp measurements.cpp blas.cpp io.cpp iram.cpp)

add_library(schwinger2d STATIC ${SCHWINGER_2D_OBJS})
target_include_directories(schwinger2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)

add_dependencies(schwinger2d Eigen)
include_directories(${EIGEN_INCLUDE_DIRS})

if(ENABLE_ALG_REMEZ)
  add_dependencies(schwinger2d algremez)
  include_directories(${ALG_REMEZ_INCLUDE_DIRS})
endif()

if(ENABLE_HDF5)
  target_link_libraries(schwinger2d PUBLIC ${HDF5_CXX_LIBRARIES})
  target_include_directories(schwinger2d PUBLIC ${HDF5_INCLUDE_DIRS})
  target_compile_definitions(schwinger2d PRIVATE ${HAVE_HDF5})
endif()

if(ENABLE_OPENMP)
  target_link_libraries(schwinger2d PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS schwinger2d)
