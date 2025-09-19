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
    core/src/zx_presets.cpp
    core/src/zx_authoring.cpp
    core/src/zx_lod.cpp
    core/src/zx_hbp.cpp
    core/src/zx_pcg.cpp
    core/src/zx_mixture.cpp
    core/src/zx_mg.cpp
    core/src/zx_integration.cpp
    core/src/zx_determinism.cpp
    core/src/zx_telemetry.cpp
    core/src/zx_checksum.cpp
    core/src/zx_residency.cpp
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
    core/src/zx_presets.cpp
    core/src/zx_authoring.cpp
    core/src/zx_lod.cpp
    core/src/zx_hbp.cpp
    core/src/zx_pcg.cpp
    core/src/zx_mixture.cpp
    core/src/zx_mg.cpp
    core/src/zx_integration.cpp
    core/src/zx_determinism.cpp
    core/src/zx_telemetry.cpp
    core/src/zx_checksum.cpp
    core/src/zx_residency.cpp
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
  add_executable(abi_conformance_test
    tests/abi_conformance_test.cpp
  )
  target_include_directories(abi_conformance_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(abi_conformance_test PRIVATE zx)
  add_test(NAME abi_conformance_test COMMAND abi_conformance_test)
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

  add_executable(anisotropy_scene_test
    tests/anisotropy_scene_test.cpp
  )
  target_include_directories(anisotropy_scene_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(anisotropy_scene_test PRIVATE zx)
  add_test(NAME anisotropy_scene_test COMMAND anisotropy_scene_test)

  add_executable(presets_test
    tests/presets_test.cpp
  )
  target_include_directories(presets_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(presets_test PRIVATE zx)
  add_test(NAME presets_test COMMAND presets_test)

  add_executable(authoring_map_test
    tests/authoring_map_test.cpp
  )
  target_include_directories(authoring_map_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(authoring_map_test PRIVATE zx)
  add_test(NAME authoring_map_test COMMAND authoring_map_test)

  add_executable(lod_filters_test
    tests/lod_filters_test.cpp
  )
  target_include_directories(lod_filters_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_filters_test PRIVATE zx)
  add_test(NAME lod_filters_test COMMAND lod_filters_test)

  add_executable(lod_alignment_test
    tests/lod_alignment_test.cpp
  )
  target_include_directories(lod_alignment_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_alignment_test PRIVATE zx)
  add_test(NAME lod_alignment_test COMMAND lod_alignment_test)

  add_executable(lod_crack_free_test
    tests/lod_crack_free_test.cpp
  )
  target_include_directories(lod_crack_free_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_crack_free_test PRIVATE zx)
  add_test(NAME lod_crack_free_test COMMAND lod_crack_free_test)

  add_executable(lod_mip_random_test
    tests/lod_mip_random_test.cpp
  )
  target_include_directories(lod_mip_random_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_mip_random_test PRIVATE zx)
  add_test(NAME lod_mip_random_test COMMAND lod_mip_random_test)

  add_executable(lod_mip_chain_test
    tests/lod_mip_chain_test.cpp
  )
  target_include_directories(lod_mip_chain_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_mip_chain_test PRIVATE zx)
  add_test(NAME lod_mip_chain_test COMMAND lod_mip_chain_test)

  add_executable(mixture_assembly_test
    tests/mixture_assembly_test.cpp
  )
  target_include_directories(mixture_assembly_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(mixture_assembly_test PRIVATE zx)
  add_test(NAME mixture_assembly_test COMMAND mixture_assembly_test)

  add_executable(bog_memory_knob_test
    tests/bog_memory_knob_test.cpp
  )
  target_include_directories(bog_memory_knob_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(bog_memory_knob_test PRIVATE zx)
  add_test(NAME bog_memory_knob_test COMMAND bog_memory_knob_test)

  add_executable(pcg_smoke_test
    tests/pcg_smoke_test.cpp
  )
  target_include_directories(pcg_smoke_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(pcg_smoke_test PRIVATE zx)
  add_test(NAME pcg_smoke_test COMMAND pcg_smoke_test)

  add_executable(mg_pcg_test
    tests/mg_pcg_test.cpp
  )
  target_include_directories(mg_pcg_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(mg_pcg_test PRIVATE zx)
  add_test(NAME mg_pcg_test COMMAND mg_pcg_test)

  add_executable(mg_prec_test
    tests/mg_prec_test.cpp
  )
  target_include_directories(mg_prec_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(mg_prec_test PRIVATE zx)
  add_test(NAME mg_prec_test COMMAND mg_prec_test)

  add_executable(dambreak_integration_test
    tests/dambreak_integration_test.cpp
  )
  target_include_directories(dambreak_integration_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(dambreak_integration_test PRIVATE zx)
  add_test(NAME dambreak_integration_test COMMAND dambreak_integration_test)

  add_executable(bogging_puddle_integration_test
    tests/bogging_puddle_integration_test.cpp
  )
  target_include_directories(bogging_puddle_integration_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(bogging_puddle_integration_test PRIVATE zx)
  add_test(NAME bogging_puddle_integration_test COMMAND bogging_puddle_integration_test)

  add_executable(lod_fallback_threshold_test
    tests/lod_fallback_threshold_test.cpp
  )
  target_include_directories(lod_fallback_threshold_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_fallback_threshold_test PRIVATE zx)
  add_test(NAME lod_fallback_threshold_test COMMAND lod_fallback_threshold_test)

  add_executable(lod_fallback_blend_test
    tests/lod_fallback_blend_test.cpp
  )
  target_include_directories(lod_fallback_blend_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_fallback_blend_test PRIVATE zx)
  add_test(NAME lod_fallback_blend_test COMMAND lod_fallback_blend_test)

  add_executable(lod_fallback_sweep_test
    tests/lod_fallback_sweep_test.cpp
  )
  target_include_directories(lod_fallback_sweep_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(lod_fallback_sweep_test PRIVATE zx)
  add_test(NAME lod_fallback_sweep_test COMMAND lod_fallback_sweep_test)

  add_executable(determinism_tests
    tests/determinism_tests.cpp
  )
  target_include_directories(determinism_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(determinism_tests PRIVATE zx)
  add_test(NAME determinism_tests COMMAND determinism_tests)

  add_executable(determinism_strict_test
    tests/determinism_strict_test.cpp
  )
  target_include_directories(determinism_strict_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(determinism_strict_test PRIVATE zx)
  add_test(NAME determinism_strict_test COMMAND determinism_strict_test)

  # Wrap the determinism script as a CTest. On Windows, use a working directory to avoid path translation issues.
  add_test(NAME determinism_script COMMAND bash wsl-determinism.sh --scene dambreak --seed 42 --fallback off)
  set_tests_properties(determinism_script PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

  add_executable(checksum_parity_test
    tests/checksum_parity_test.cpp
  )
  target_include_directories(checksum_parity_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(checksum_parity_test PRIVATE zx)
  add_test(NAME checksum_parity_test COMMAND checksum_parity_test)

  add_executable(determinism_ordering_test
    tests/determinism_ordering_test.cpp
  )
  target_include_directories(determinism_ordering_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(determinism_ordering_test PRIVATE zx)
  add_test(NAME determinism_ordering_test COMMAND determinism_ordering_test)

  add_executable(residency_hysteresis_test
    tests/residency_hysteresis_test.cpp
  )
  target_include_directories(residency_hysteresis_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(residency_hysteresis_test PRIVATE zx)
  add_test(NAME residency_hysteresis_test COMMAND residency_hysteresis_test)

  add_executable(residency_prefetch_test
    tests/residency_prefetch_test.cpp
  )
  target_include_directories(residency_prefetch_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(residency_prefetch_test PRIVATE zx)
  add_test(NAME residency_prefetch_test COMMAND residency_prefetch_test)

  add_executable(residency_pinning_test
    tests/residency_pinning_test.cpp
  )
  target_include_directories(residency_pinning_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(residency_pinning_test PRIVATE zx)
  add_test(NAME residency_pinning_test COMMAND residency_pinning_test)

  add_executable(residency_churn_test
    tests/residency_churn_test.cpp
  )
  target_include_directories(residency_churn_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(residency_churn_test PRIVATE zx)
  add_test(NAME residency_churn_test COMMAND residency_churn_test)

  add_executable(hbp_mu_eff_test
    tests/hbp_mu_eff_test.cpp
  )
  target_include_directories(hbp_mu_eff_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(hbp_mu_eff_test PRIVATE zx)
  add_test(NAME hbp_mu_eff_test COMMAND hbp_mu_eff_test)

  add_executable(hbp_viscous_update_test
    tests/hbp_viscous_update_test.cpp
  )
  target_include_directories(hbp_viscous_update_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(hbp_viscous_update_test PRIVATE zx)
  add_test(NAME hbp_viscous_update_test COMMAND hbp_viscous_update_test)

  add_executable(hbp_clamp_policy_test
    tests/hbp_clamp_policy_test.cpp
  )
  target_include_directories(hbp_clamp_policy_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(hbp_clamp_policy_test PRIVATE zx)
  add_test(NAME hbp_clamp_policy_test COMMAND hbp_clamp_policy_test)

  add_executable(hbp_stability_bound_test
    tests/hbp_stability_bound_test.cpp
  )
  target_include_directories(hbp_stability_bound_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(hbp_stability_bound_test PRIVATE zx)
  add_test(NAME hbp_stability_bound_test COMMAND hbp_stability_bound_test)

  add_executable(telemetry_schema_test
    tests/telemetry_schema_test.cpp
  )
  target_include_directories(telemetry_schema_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(telemetry_schema_test PRIVATE zx)
  add_test(NAME telemetry_schema_test COMMAND telemetry_schema_test)

  add_executable(telemetry_depth_test
    tests/telemetry_depth_test.cpp
  )
  target_include_directories(telemetry_depth_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core/include)
  target_link_libraries(telemetry_depth_test PRIVATE zx)
  add_test(NAME telemetry_depth_test COMMAND telemetry_depth_test)
endif()


