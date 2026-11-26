# CMake for the rssn C++ client (header-only).
#
# The client is two headers plus the generated operator facade; it links against
# the rssn C library, which the repository builds with
# `bash scripts/build_ffi.sh`.
#
#   cmake -S bindings/cpp -B build/cpp \
#         -DRSSN_ROOT=<repo> \
#         -DRSSN_LIBRARY=<repo>/target/ffi/librssn.a
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

add_library(rssn_client INTERFACE)
target_include_directories(rssn_client INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${RSSN_ROOT}"
)
target_compile_features(rssn_client INTERFACE cxx_std_17)

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

if(DEFINED RSSN_LIBRARY)
    target_link_libraries(rssn_smoke PRIVATE "${RSSN_LIBRARY}")
    target_link_libraries(rssn_conformance PRIVATE "${RSSN_LIBRARY}")
    target_link_libraries(rssn_callbacks PRIVATE "${RSSN_LIBRARY}")
    target_link_libraries(rssn_workflow PRIVATE "${RSSN_LIBRARY}")

    # A shared library must be findable at run time: put it next to the test
    # binaries on Windows, and on the rpath elsewhere.
    get_filename_component(RSSN_LIBRARY_DIR "${RSSN_LIBRARY}" DIRECTORY)
    if(WIN32 AND EXISTS "${RSSN_LIBRARY_DIR}/rssn.dll")
        foreach(target rssn_smoke rssn_conformance rssn_callbacks rssn_workflow)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        "${RSSN_LIBRARY_DIR}/rssn.dll" "$<TARGET_FILE_DIR:${target}>")
        endforeach()
    elseif(UNIX)
        set_target_properties(rssn_smoke rssn_conformance rssn_callbacks rssn_workflow PROPERTIES
            BUILD_RPATH "${RSSN_LIBRARY_DIR}"
            INSTALL_RPATH "${RSSN_LIBRARY_DIR}")
    endif()
else()
    message(WARNING "RSSN_LIBRARY is not set: the targets compile but will not link")
endif()

enable_testing()
add_test(NAME smoke COMMAND rssn_smoke)
add_test(NAME conformance
         COMMAND rssn_conformance "${CMAKE_CURRENT_SOURCE_DIR}/../tests/cases.tsv")
add_test(NAME workflow COMMAND rssn_workflow)
