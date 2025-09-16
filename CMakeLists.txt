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
    core/src/zx_contact_ref.cpp
    core/src/zx_writeback.cpp
    core/src/zx_constitutive_ref.cpp
    core/src/zx_validation.cpp
    core/src/zx_heave_validation.cpp
  )
else()
  add_library(zx STATIC
    core/src/zx_api.cpp
    core/src/zx_tiles.cpp
    core/src/zx_apic_ref.cpp
    core/src/zx_contact_ref.cpp
    core/src/zx_writeback.cpp
    core/src/zx_constitutive_ref.cpp
    core/src/zx_validation.cpp
    core/src/zx_heave_validation.cpp
  )
endif()

add_executable(zx_cli tools/zx_cli.cpp)
target_include_directories(zx_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
target_link_libraries(zx_cli PRIVATE zx)

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

  add_executable(constitutive_ref_test
    tests/constitutive_ref_test.cpp
  )
  target_include_directories(constitutive_ref_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(constitutive_ref_test PRIVATE zx)
  add_test(NAME constitutive_ref_test COMMAND constitutive_ref_test)

  add_executable(constitutive_plastic_work_test
    tests/constitutive_plastic_work_test.cpp
  )
  target_include_directories(constitutive_plastic_work_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(constitutive_plastic_work_test PRIVATE zx)
  add_test(NAME constitutive_plastic_work_test COMMAND constitutive_plastic_work_test)

  add_executable(mc_return_map_test
    tests/mc_return_map_test.cpp
  )
  target_include_directories(mc_return_map_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(mc_return_map_test PRIVATE zx)
  add_test(NAME mc_return_map_test COMMAND mc_return_map_test)

  add_executable(validation_tests
    tests/validation_tests.cpp
  )
  target_include_directories(validation_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(validation_tests PRIVATE zx)
  add_test(NAME validation_tests COMMAND validation_tests)

  add_executable(norsand_test
    tests/norsand_test.cpp
  )
  target_include_directories(norsand_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(norsand_test PRIVATE zx)
  add_test(NAME norsand_test COMMAND norsand_test)

  add_executable(heave_validation_test
    tests/heave_validation_test.cpp
  )
  target_include_directories(heave_validation_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(heave_validation_test PRIVATE zx)
  add_test(NAME heave_validation_test COMMAND heave_validation_test)

  add_executable(anisotropic_contact_test
    tests/anisotropic_contact_test.cpp
  )
  target_include_directories(anisotropic_contact_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(anisotropic_contact_test PRIVATE zx)
  add_test(NAME anisotropic_contact_test COMMAND anisotropic_contact_test)

  add_executable(wheel_heave_test
    tests/wheel_heave_test.cpp
  )
  target_include_directories(wheel_heave_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(wheel_heave_test PRIVATE zx)
  add_test(NAME wheel_heave_test COMMAND wheel_heave_test)

  add_executable(foot_heave_test
    tests/foot_heave_test.cpp
  )
  target_include_directories(foot_heave_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(foot_heave_test PRIVATE zx)
  add_test(NAME foot_heave_test COMMAND foot_heave_test)
endif()


