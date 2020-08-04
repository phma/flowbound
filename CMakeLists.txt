project(flowbound)
cmake_minimum_required(VERSION 3.8.0)
add_executable(flowbound flowbound.cpp eisenstein.cpp flowsnake.cpp)

include(TestBigEndian)
include(CheckTypeSize)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17) # appeared in CMake 3.8
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(SHARE_DIR ${CMAKE_INSTALL_PREFIX}/share/flowbound)

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." FALSE)
if (${FORCE_COLORED_OUTPUT})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options (-fdiagnostics-color=always)
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options (-fcolor-diagnostics)
  endif ()
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(GMPXX REQUIRED)
target_link_libraries(flowbound ${GMP_LIBRARY} ${GMPXX_LIBRARIES})

install(TARGETS flowbound DESTINATION bin)

test_big_endian(BIGENDIAN)
check_type_size("int" INT)
check_type_size("int *" INT_POINTER)

set(FLOWBOUND_MAJOR_VERSION 0)
set(FLOWBOUND_MINOR_VERSION 1)
set(FLOWBOUND_PATCH_VERSION 0)
set(FLOWBOUND_VERSION ${FLOWBOUND_MAJOR_VERSION}.${FLOWBOUND_MINOR_VERSION}.${FLOWBOUND_PATCH_VERSION})

set(CPACK_PACKAGE_VERSION_MAJOR ${FLOWBOUND_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${FLOWBOUND_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${FLOWBOUND_PATCH_VERSION})
set(CPACK_SOURCE_IGNORE_FILES /\\\\.git;.*~)
include(CPack)

include(CTest)
