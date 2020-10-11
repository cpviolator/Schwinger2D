cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the project name and version
project(Schwinger2D VERSION 1.0.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} show
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE IS_GIT_REPOSIITORY
    OUTPUT_QUIET ERROR_QUIET)
  if(${IS_GIT_REPOSIITORY} EQUAL 0)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --abbrev=0
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GITTAG OUTPUT_STRIP_TRAILING_WHITESPACE)
    # we use git rev-list and pipe that through wc here. Newer git versions support --count as option to rev-list but
    # that might not always be available
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-list ${GITTAG}..HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      COMMAND wc -l
      OUTPUT_VARIABLE GITCOUNT OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --long --dirty
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GITVERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
endif(GIT_FOUND)


add_subdirectory(lib)
add_subdirectory(wilson)
