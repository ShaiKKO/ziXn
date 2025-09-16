##!
# \file CMakeLists.txt
# \brief Top-level CMake build for ziXn.
# \author Colin Macritchie (Ripple Group, LLC)
# \version 1.0.0
# \date 2025-09-16
# \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
##

cmake_minimum_required(VERSION 3.20)
project(ziXn LANGUAGES C CXX VERSION 1.0.0)

option(ZX_BUILD_SHARED "Build ziXn as a shared library" ON)
option(ZX_ENABLE_TESTS "Build tests" ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ZX_BUILD_SHARED)
  add_library(zx SHARED
    core/src/zx_api.cpp
    core/src/zx_tiles.cpp
    core/src/zx_apic_ref.cpp
    core/src/zx_writeback.cpp
  )
else()
  add_library(zx STATIC
    core/src/zx_api.cpp
    core/src/zx_tiles.cpp
    core/src/zx_apic_ref.cpp
    core/src/zx_writeback.cpp
  )
endif()

target_include_directories(zx
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/core/include
)

target_compile_definitions(zx
  PRIVATE ZX_BUILD
)

if(MSVC)
  target_compile_options(zx PRIVATE /W4 /permissive-)
else()
  target_compile_options(zx PRIVATE -Wall -Wextra -Wpedantic)
endif()

set_target_properties(zx PROPERTIES
  OUTPUT_NAME zx
  PREFIX ""
)

install(TARGETS zx
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(DIRECTORY core/include/ DESTINATION include)

if(ZX_ENABLE_TESTS)
  include(CTest)
  enable_testing()
  add_executable(zx_abi_smoke
    tests/abi_smoke.cpp
  )
  target_include_directories(zx_abi_smoke PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(zx_abi_smoke PRIVATE zx)
  add_test(NAME zx_abi_smoke COMMAND zx_abi_smoke)

  add_executable(apic_basis_test
    tests/apic_basis_test.cpp
  )
  add_test(NAME apic_basis_test COMMAND apic_basis_test)

  add_executable(apic_ref_test
    tests/apic_ref_test.cpp
  )
  target_include_directories(apic_ref_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(apic_ref_test PRIVATE zx)
  add_test(NAME apic_ref_test COMMAND apic_ref_test)

  add_executable(counters_smoke
    tests/counters_smoke.cpp
  )
  target_include_directories(counters_smoke PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(counters_smoke PRIVATE zx)
  add_test(NAME counters_smoke COMMAND counters_smoke)

  add_executable(writeback_smoke
    tests/writeback_smoke.cpp
  )
  target_include_directories(writeback_smoke PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(writeback_smoke PRIVATE zx)
  add_test(NAME writeback_smoke COMMAND writeback_smoke)
endif()


