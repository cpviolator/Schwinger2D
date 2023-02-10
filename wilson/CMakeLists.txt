add_executable(wilson2D wilson2D.cpp)
target_include_directories(wilson2D PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(wilson2D schwinger2d)

include_directories(${EIGEN_INCLUDE_DIRS})

if(ENABLE_ALG_REMEZ)
  target_link_libraries(wilson2D algremez ${GMP_LIB} ${MPFR_LIBRARIES})
  include_directories(${ALG_REMEZ_INCLUDE_DIRS})
endif()

if(ENABLE_FEAST)
  #target_link_libraries(wilson2D feast ${MKL_LIB_DIRS}/libmkl_gf_lp64.so.2 ${MKL_LIB_DIRS}/libmkl_gnu_thread.so.2 ${MKL_LIB_DIRS}/libmkl_core.so.2 m gfortran)
  #target_link_libraries(wilson2D feast m gfortran)
  #include_directories(${FEAST_INCLUDE_DIRS})
endif()


