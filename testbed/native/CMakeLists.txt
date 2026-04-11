cmake_minimum_required(VERSION 3.21)
project(testbed_native C CXX)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use dynamic CRT everywhere (required for DLL builds, must match Jolt)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
set(USE_STATIC_MSVC_RUNTIME_LIBRARY OFF CACHE BOOL "" FORCE)

# Output DLLs to bin/ for the .NET app to find
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# nudge shared library
# ---------------------------------------------------------------------------

add_library(nudge_shared SHARED nudge_dll.c)
target_include_directories(nudge_shared PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
set_target_properties(nudge_shared PROPERTIES OUTPUT_NAME "nudge")

if(MSVC)
	target_compile_options(nudge_shared PRIVATE /experimental:c11atomics /Ob3)
endif()

# ---------------------------------------------------------------------------
# JoltPhysics (fetched from source, wrapped with flat C exports)
# ---------------------------------------------------------------------------

include(FetchContent)
FetchContent_Declare(
	JoltPhysics
	GIT_REPOSITORY https://github.com/jrouwe/JoltPhysics.git
	GIT_TAG master
	GIT_SHALLOW TRUE
	SOURCE_SUBDIR Build
)
FetchContent_MakeAvailable(JoltPhysics)

add_library(jolt_shared SHARED jolt_dll.cpp)
target_link_libraries(jolt_shared PRIVATE Jolt::Jolt)
set_target_properties(jolt_shared PROPERTIES OUTPUT_NAME "jolt")
