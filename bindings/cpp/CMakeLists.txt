# CMake for the rssn C++ client (header-only).
#
# The client is two headers plus the generated operator facade; it links against
# the rssn C library, which the repository builds with
# `bazel build //:ffi_distribution`.
#
#   cmake -S bindings/cpp -B build/cpp \
#         -DRSSN_ROOT=<repo> \
#         -DRSSN_LIBRARY=<repo>/bazel-bin/ffi_distribution/ffi/librssn.a
#   cmake --build build/cpp
#   ctest --test-dir build/cpp --output-on-failure
#
# `RSSN_LIBRARY` is the import library (`rssn.dll.lib` on MSVC,
# `librssn.dll.a` on MinGW) or the static archive, whichever matches the
# toolchain. On MSVC an import library is produced next to `rssn.dll`; on MinGW
# generate one from the DLL with
#
#   gendef rssn.dll && dlltool -d rssn.def -l librssn.dll.a
#
cmake_minimum_required(VERSION 3.16)
project(rssn_cpp_client LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED RSSN_ROOT)
    set(RSSN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
endif()
if(NOT EXISTS "${RSSN_ROOT}/rssn.h")
    message(FATAL_ERROR "RSSN_ROOT must be the rssn repository root (it holds rssn.h); got ${RSSN_ROOT}")
endif()

# Use the same shared-library file or directory as the dynamic bindings.
# An explicit RSSN_LIBRARY still selects the file used at link time.
if(NOT RSSN_LIBRARY)
    if(WIN32)
        if(MSVC)
            set(rssn_link_names rssn.dll.lib rssn.lib)
        else()
            set(rssn_link_names librssn.dll.a librssn.a)
        endif()
    elseif(APPLE)
        set(rssn_link_names librssn.dylib librssn.a)
    else()
        set(rssn_link_names librssn.so librssn.so.0 librssn.a)
    endif()
    if(NOT "$ENV{RSSN_LIB}" STREQUAL "")
        get_filename_component(rssn_location "$ENV{RSSN_LIB}" ABSOLUTE BASE_DIR "${RSSN_ROOT}")
        if(IS_DIRECTORY "${rssn_location}")
            set(rssn_search "${rssn_location}")
        elseif(WIN32 AND EXISTS "${rssn_location}")
            get_filename_component(rssn_search "${rssn_location}" DIRECTORY)
            set(rssn_dll "${rssn_location}")
        elseif(EXISTS "${rssn_location}")
            set(RSSN_LIBRARY "${rssn_location}")
        endif()
    else()
        list(APPEND rssn_search "${RSSN_ROOT}")
        if(EXISTS "${RSSN_ROOT}/MODULE.bazel")
            foreach(target ffi_distribution ffi ffi_distribution_portable ffi_dev ffi_distribution_dev ffi_static ffi_static_dev)
                list(APPEND rssn_search "${RSSN_ROOT}/bazel-bin/${target}/ffi")
            endforeach()
        endif()
        foreach(profile debug ffi release)
            foreach(directory target target-bindings target-ffi dist/rssn)
                list(APPEND rssn_search "${RSSN_ROOT}/${directory}/${profile}")
            endforeach()
        endforeach()
    endif()
    foreach(directory IN LISTS rssn_search)
        foreach(name IN LISTS rssn_link_names)
            if(EXISTS "${directory}/${name}" AND NOT IS_DIRECTORY "${directory}/${name}")
                set(RSSN_LIBRARY "${directory}/${name}")
                break()
            endif()
        endforeach()
        if(RSSN_LIBRARY)
            break()
        endif()
    endforeach()
endif()
if(RSSN_LIBRARY)
    get_filename_component(RSSN_LIBRARY "${RSSN_LIBRARY}" ABSOLUTE BASE_DIR "${RSSN_ROOT}")
endif()
if(NOT RSSN_LIBRARY OR NOT EXISTS "${RSSN_LIBRARY}" OR IS_DIRECTORY "${RSSN_LIBRARY}")
    message(FATAL_ERROR
        "RSSN library for this C++ toolchain not found. Build it with bazel build //:ffi_dev, "
        "set RSSN_LIB to the shared library or its directory, or pass -DRSSN_LIBRARY=<link library>.")
endif()
message(STATUS "RSSN link library: ${RSSN_LIBRARY}")

add_library(rssn_client INTERFACE)
target_include_directories(rssn_client INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${RSSN_ROOT}"
)
target_compile_features(rssn_client INTERFACE cxx_std_17)
get_filename_component(rssn_library_extension "${RSSN_LIBRARY}" LAST_EXT)
string(TOLOWER "${rssn_library_extension}" rssn_library_extension)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND rssn_library_extension STREQUAL ".a")
    find_package(Threads REQUIRED)
    target_link_libraries(rssn_client INTERFACE "${RSSN_LIBRARY}" m Threads::Threads ${CMAKE_DL_LIBS})
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND RSSN_LIBRARY MATCHES "[.]so([.]0)?$")
    add_library(rssn_native SHARED IMPORTED)
    set_target_properties(rssn_native PROPERTIES
        IMPORTED_LOCATION "${RSSN_LIBRARY}"
        IMPORTED_NO_SONAME TRUE)
    target_link_libraries(rssn_client INTERFACE rssn_native)
else()
    target_link_libraries(rssn_client INTERFACE "${RSSN_LIBRARY}")
endif()

option(RSSN_CPP_BUILD_TESTS "Build the RSSN C++ SDK tests and examples" ON)
if(NOT RSSN_CPP_BUILD_TESTS)
    return()
endif()

add_executable(rssn_smoke tests/smoke.cpp)
target_link_libraries(rssn_smoke PRIVATE rssn_client)

add_executable(rssn_conformance "${CMAKE_CURRENT_SOURCE_DIR}/../tests/conformance.cpp")
target_link_libraries(rssn_conformance PRIVATE rssn_client)

# The workflows `README.md` shows, as a program (`bindings/tests/run_readme_examples.py`).
add_executable(rssn_workflow tests/workflow.cpp)
target_link_libraries(rssn_workflow PRIVATE rssn_client)

# The host-callback example (`bindings/examples/cpp/callbacks.cpp`).
add_executable(rssn_callbacks "${CMAKE_CURRENT_SOURCE_DIR}/../examples/cpp/callbacks.cpp")
target_link_libraries(rssn_callbacks PRIVATE rssn_client)

add_executable(rssn_exact_tensor "${CMAKE_CURRENT_SOURCE_DIR}/../examples/cpp/exact_tensor.cpp")
target_link_libraries(rssn_exact_tensor PRIVATE rssn_client)

# A shared library must be findable at run time for the test executables.
get_filename_component(RSSN_LIBRARY_DIR "${RSSN_LIBRARY}" DIRECTORY)
if(WIN32)
    if(NOT rssn_dll)
        set(rssn_dll "${RSSN_LIBRARY_DIR}/rssn.dll")
    endif()
    if(EXISTS "${rssn_dll}")
        foreach(target rssn_smoke rssn_conformance rssn_callbacks rssn_workflow rssn_exact_tensor)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        "${rssn_dll}" "$<TARGET_FILE_DIR:${target}>"
                VERBATIM)
        endforeach()
    endif()
elseif(UNIX AND NOT rssn_library_extension STREQUAL ".a")
    set_target_properties(rssn_smoke rssn_conformance rssn_callbacks rssn_workflow rssn_exact_tensor PROPERTIES
        BUILD_RPATH "${RSSN_LIBRARY_DIR}")
endif()

enable_testing()
add_test(NAME smoke COMMAND rssn_smoke)
add_test(NAME callbacks COMMAND rssn_callbacks)
add_test(NAME conformance
         COMMAND rssn_conformance "${CMAKE_CURRENT_SOURCE_DIR}/../tests/cases.tsv")
add_test(NAME workflow COMMAND rssn_workflow)
add_test(NAME exact_tensor COMMAND rssn_exact_tensor)
