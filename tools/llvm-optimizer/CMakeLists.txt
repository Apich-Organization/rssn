cmake_minimum_required(VERSION 3.20)
project(RssnLLVMOptimizer LANGUAGES CXX)

set(RSSN_RUSTC "rustc" CACHE FILEPATH "rustc that will load the plugin")
set(RSSN_LLVM_CONFIG "llvm-config-23" CACHE FILEPATH "LLVM SDK llvm-config")
option(RSSN_TEST_RUSTC "Test plugin loading in native rustc" ON)

execute_process(COMMAND "${RSSN_RUSTC}" -vV
  OUTPUT_VARIABLE rust_version COMMAND_ERROR_IS_FATAL ANY)
string(REGEX MATCH "LLVM version: (23\\.[0-9]+\\.[0-9]+)" rust_llvm "${rust_version}")
set(rust_llvm_version "${CMAKE_MATCH_1}")
if(NOT rust_llvm_version)
  message(FATAL_ERROR "This plugin requires a rustc built with LLVM 23")
endif()

execute_process(COMMAND "${RSSN_LLVM_CONFIG}" --version
  OUTPUT_VARIABLE sdk_version OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
if(NOT sdk_version STREQUAL rust_llvm_version)
  message(FATAL_ERROR "rustc uses LLVM ${rust_llvm_version}, SDK is ${sdk_version}")
endif()
add_library(RssnLLVMOptimizer MODULE Pressure.cpp SpillPressure.cpp Motion.cpp LoadFusion.cpp ScalarLoadPack.cpp SwitchLoadReuse.cpp MemoryModel.cpp CallABI.cpp LibmPressure.cpp LibmReuse.cpp LibmTailCutoff.cpp LibmUnitLog.cpp PowExponentPack.cpp LoopExitCompression.cpp LoopPreheaderPressure.cpp FPClassInterval.cpp FPEndpointGuards.cpp FPRegisterBank.cpp PowUnitBase.cpp)
target_sources(RssnLLVMOptimizer PRIVATE EdgeRestore.cpp MachineStackTraffic.cpp)
# rustc's LLVM builds without RTTI. A distro SDK can otherwise emit imports
# for the machine pass and diagnostic base-class typeinfo.
if(MSVC)
  set_source_files_properties(EdgeRestore.cpp MachineStackTraffic.cpp PROPERTIES COMPILE_OPTIONS /GR-)
else()
  set_source_files_properties(EdgeRestore.cpp MachineStackTraffic.cpp PROPERTIES COMPILE_OPTIONS -fno-rtti)
endif()
target_compile_features(RssnLLVMOptimizer PRIVATE cxx_std_17)
set_target_properties(RssnLLVMOptimizer PROPERTIES PREFIX "" CXX_VISIBILITY_PRESET hidden)
if(MSVC)
  execute_process(COMMAND "${RSSN_LLVM_CONFIG}" --includedir
    OUTPUT_VARIABLE sdk_include OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
  execute_process(COMMAND "${RSSN_LLVM_CONFIG}" --has-rtti
    OUTPUT_VARIABLE sdk_rtti OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
  target_include_directories(RssnLLVMOptimizer SYSTEM PRIVATE "${sdk_include}")
  if(sdk_rtti STREQUAL "NO")
    target_compile_options(RssnLLVMOptimizer PRIVATE /GR-)
  endif()
else()
  execute_process(COMMAND "${RSSN_LLVM_CONFIG}" --cxxflags
    OUTPUT_VARIABLE sdk_flags OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
  separate_arguments(sdk_flags NATIVE_COMMAND "${sdk_flags}")
  target_compile_options(RssnLLVMOptimizer PRIVATE ${sdk_flags})
endif()

# Use the host process's LLVM so pass and analysis identities stay shared.
if(WIN32)
  set(RSSN_LLVM_IMPORT_LIBRARY "" CACHE FILEPATH "Import library for the host's LLVM C++ exports")
  if(NOT EXISTS "${RSSN_LLVM_IMPORT_LIBRARY}")
    message(FATAL_ERROR "Set RSSN_LLVM_IMPORT_LIBRARY to the import library for the executable or DLL exporting the host's LLVM C++ API")
  endif()
  target_link_libraries(RssnLLVMOptimizer PRIVATE "${RSSN_LLVM_IMPORT_LIBRARY}")
elseif(APPLE)
  target_link_options(RssnLLVMOptimizer PRIVATE "LINKER:-undefined,dynamic_lookup")
endif()

include(CTest)
if(BUILD_TESTING AND NOT CMAKE_CROSSCOMPILING)
  # Retain CTest names for existing callers while Bazel supplies the sources,
  # plugin and test arguments.
  find_program(bazel NAMES bazel bazelisk)
  if(NOT bazel)
    set(bazel bazel)
  endif()
  find_program(sdk_program NAMES "${RSSN_LLVM_CONFIG}" REQUIRED)
  find_program(rustc_program NAMES "${RSSN_RUSTC}" REQUIRED)
  set(bazel_flags "--repo_env=RSSN_LLVM_CONFIG=${sdk_program}"
                  "--repo_env=RUSTC=${rustc_program}")
  if(WIN32)
    list(APPEND bazel_flags "--repo_env=RSSN_LLVM_IMPORT_LIBRARY=${RSSN_LLVM_IMPORT_LIBRARY}")
  endif()
  add_test(NAME llvm-optimizer-ir
    COMMAND "${bazel}" test ${bazel_flags} //:llvm_optimizer_tools_test
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../..")
  if(RSSN_TEST_RUSTC)
    add_test(NAME llvm-optimizer-rustc
      COMMAND "${bazel}" test ${bazel_flags} //:llvm_optimizer_rustc_test
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../..")
  endif()
endif()
