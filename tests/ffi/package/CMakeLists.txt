cmake_minimum_required(VERSION 3.20)
project(rssn_packaged_clients LANGUAGES C CXX)

if(RSSN_LIBRARY_DIR)
    if(RSSN_ARCHIVE OR NOT RSSN_LINKAGE MATCHES "^(shared|static)$")
        message(FATAL_ERROR "Select a library directory and a shared or static linkage")
    endif()
    if(RSSN_LINKAGE STREQUAL "static" AND NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message(FATAL_ERROR "The static client test is configured for Linux")
    endif()
    get_filename_component(package "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
    set(libraries "${RSSN_LIBRARY_DIR}")
else()
    if(NOT RSSN_ARCHIVE OR NOT EXISTS "${RSSN_ARCHIVE}" OR NOT EXISTS "${RSSN_ARCHIVE}.sha256")
        message(FATAL_ERROR "RSSN_ARCHIVE must point to a release archive with its .sha256 companion")
    endif()
    get_filename_component(archive_name "${RSSN_ARCHIVE}" NAME)
    string(REGEX REPLACE "\\.tar\\.gz$" "" package_name "${archive_name}")
    if(archive_name STREQUAL package_name)
        message(FATAL_ERROR "The release archive must use the .tar.gz format")
    endif()
    file(SHA256 "${RSSN_ARCHIVE}" actual_hash)
    file(READ "${RSSN_ARCHIVE}.sha256" checksum)
    string(STRIP "${checksum}" checksum)
    string(REGEX MATCH "^([0-9A-Fa-f]+)[ \t]+([^\r\n]+)$" checksum_line "${checksum}")
    string(TOLOWER "${CMAKE_MATCH_1}" expected_hash)
    if(NOT expected_hash STREQUAL actual_hash OR NOT CMAKE_MATCH_2 STREQUAL archive_name)
        message(FATAL_ERROR "The packaged archive does not match its named SHA-256 digest")
    endif()
    file(ARCHIVE_EXTRACT INPUT "${RSSN_ARCHIVE}" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/archive")
    set(package "${CMAKE_CURRENT_BINARY_DIR}/archive/${package_name}")
    set(libraries "${package}")
    foreach(sdk_entry
            bindings/cpp/CMakeLists.txt
            bindings/cpp/include/rssn/client.hpp
            bindings/julia/Project.toml
            bindings/julia/src/RSSN.jl
            bindings/python/pyproject.toml
            bindings/python/rssn/__init__.py
            bindings/r/DESCRIPTION
            bindings/r/R/zzz.R
            bindings/r/src/rssn_r.c)
        if(NOT EXISTS "${package}/${sdk_entry}")
            message(FATAL_ERROR "Missing packaged SDK file ${sdk_entry}")
        endif()
    endforeach()
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

add_library(rssn_shared SHARED IMPORTED)
if(WIN32)
    set_target_properties(rssn_shared PROPERTIES
        IMPORTED_LOCATION "${libraries}/rssn.dll"
        IMPORTED_IMPLIB "${libraries}/rssn.dll.lib")
elseif(APPLE)
    set_target_properties(rssn_shared PROPERTIES IMPORTED_LOCATION "${libraries}/librssn.dylib")
else()
    # A Rust cdylib can omit DT_SONAME. Search by library name so the
    # consumers keep using this package's RPATH from any working directory.
    set_target_properties(rssn_shared PROPERTIES
        IMPORTED_LOCATION "${libraries}/librssn.so"
        IMPORTED_NO_SONAME TRUE)
endif()

if(RSSN_ARCHIVE AND NOT CMAKE_CROSSCOMPILING)
    find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
    add_test(NAME packaged_python COMMAND "${Python3_EXECUTABLE}" -I -B
        "${CMAKE_CURRENT_SOURCE_DIR}/sdk_smoke.py" "${package}" "$<TARGET_FILE:rssn_shared>")
endif()

function(exercise_package kind target)
    add_executable("${kind}_smoke" "${package}/examples/c/smoke.c")
    add_executable("${kind}_exact_key" "${package}/examples/c/operator_contract.c")
    add_executable("${kind}_headers" ../headers.cpp)
    add_executable("${kind}_cpp_workflow" "${package}/bindings/cpp/tests/workflow.cpp")
    add_executable("${kind}_features" backend_features.cpp)
    target_include_directories("${kind}_cpp_workflow" PRIVATE "${package}/bindings/cpp/include")
    target_include_directories("${kind}_features" PRIVATE "${package}/bindings/cpp/include")
    foreach(client smoke exact_key headers cpp_workflow features)
        set(executable "${kind}_${client}")
        target_include_directories("${executable}" PRIVATE "${package}")
        target_link_libraries("${executable}" PRIVATE "${target}")
        if(MSVC)
            target_compile_options("${executable}" PRIVATE /utf-8)
        endif()
        if(WIN32 AND kind STREQUAL "shared")
            add_custom_command(TARGET "${executable}" POST_BUILD
                COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                    "$<TARGET_FILE:${target}>" "$<TARGET_FILE_DIR:${executable}>")
        endif()
        if(client STREQUAL "exact_key")
            add_test(NAME "${executable}" COMMAND "${CMAKE_COMMAND}"
                "-DCLIENT=$<TARGET_FILE:${executable}>" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_exact_key.cmake")
        elseif(client STREQUAL "features")
            add_test(NAME "${kind}_manifest" COMMAND "${CMAKE_COMMAND}"
                "-DCLIENT=$<TARGET_FILE:${executable}>" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_capabilities.cmake")
            add_test(NAME "${kind}_backend" COMMAND "${executable}" backend)
            add_test(NAME "${kind}_npy" COMMAND "${executable}" npy)
        else()
            add_test(NAME "${executable}" COMMAND "${executable}")
        endif()
    endforeach()
endfunction()

if(NOT RSSN_LINKAGE STREQUAL "static")
    exercise_package(shared rssn_shared)
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND NOT RSSN_LINKAGE STREQUAL "shared")
    set(RSSN_ROOT "${package}")
    set(RSSN_LIBRARY "${libraries}/librssn.a")
    set(RSSN_CPP_BUILD_TESTS OFF)
    add_subdirectory("${package}/bindings/cpp" "${CMAKE_CURRENT_BINARY_DIR}/cpp-sdk" EXCLUDE_FROM_ALL)
    exercise_package(static rssn_client)
endif()
