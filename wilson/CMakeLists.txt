add_executable(wilson2D wilson2D.cpp)
target_link_libraries(wilson2D schwinger2d algremez ${GMP_LIB} ${MPFR_LIBRARIES})
target_include_directories(wilson2D PUBLIC ${PROJECT_SOURCE_DIR}/include)
include_directories(${EIGEN_INCLUDE_DIRS} ${ALG_REMEZ_INCLUDE_DIRS})
