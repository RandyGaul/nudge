cmake_minimum_required(VERSION 3.21)
project(nudge C CXX)

# CMake's CMAKE_C_STANDARD=23 isn't recognized for MSVC, so we pass the flag
# directly. Non-MSVC compilers get C23 via CMAKE_C_STANDARD as usual.
if(NOT MSVC)
	set(CMAKE_C_STANDARD 23)
	set(CMAKE_C_STANDARD_REQUIRED ON)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The GUI app, ImGui, and the RPM viewer need SDL3 + OpenGL + (for the
# viewer) Win32. CI and headless builds just need nudge_tests, which has
# zero external deps. Toggle NUDGE_BUILD_APP=OFF to skip the heavy parts.
option(NUDGE_BUILD_APP "Build the GUI app + ImGui + remote viewer" ON)

# ----- always-on: the headless test binary ------------------------------------
add_executable(nudge_tests src/test_main.c)

# pthreads + libm on Unix-like platforms (empty on Windows, which uses Win32
# directly and links the math runtime as part of the CRT).
if(NOT WIN32)
	set(THREADS_PREFER_PTHREAD_FLAG ON)
	find_package(Threads REQUIRED)
	target_link_libraries(nudge_tests Threads::Threads m)
endif()

if(MSVC)
	target_compile_options(nudge_tests PRIVATE /std:clatest /experimental:c11atomics /Ob3)
	target_link_options(nudge_tests PRIVATE /STACK:8388608)
endif()

if(NUDGE_BUILD_APP)
	include(FetchContent)
	FetchContent_Declare(
		SDL3
		GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
		GIT_TAG main
		GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(SDL3)

	# Dear ImGui (compiled as C++ separately from the unity C build)
	set(NUDGE_IMGUI_DIR ${CMAKE_SOURCE_DIR}/libraries/imgui)
	add_library(nudge_imgui STATIC
		${NUDGE_IMGUI_DIR}/imgui.cpp
		${NUDGE_IMGUI_DIR}/imgui_draw.cpp
		${NUDGE_IMGUI_DIR}/imgui_tables.cpp
		${NUDGE_IMGUI_DIR}/imgui_widgets.cpp
		${NUDGE_IMGUI_DIR}/imgui_demo.cpp
		${NUDGE_IMGUI_DIR}/imgui_impl_sdl3.cpp
		${NUDGE_IMGUI_DIR}/imgui_impl_opengl3.cpp
		${NUDGE_IMGUI_DIR}/dcimgui.cpp
		${NUDGE_IMGUI_DIR}/dcimgui_impl_sdl3.cpp
		${NUDGE_IMGUI_DIR}/dcimgui_impl_opengl3.cpp
	)
	target_include_directories(nudge_imgui PUBLIC ${NUDGE_IMGUI_DIR})
	target_link_libraries(nudge_imgui PUBLIC SDL3::SDL3)

	set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT nudge)
	add_executable(nudge src/main.c)
	# OpenGL lib name differs per platform. On platforms we don't support the
	# GUI on, this just fails at link time (harmless if NUDGE_BUILD_APP=OFF).
	if(WIN32)
		target_link_libraries(nudge nudge_imgui opengl32)
	else()
		find_package(OpenGL)
		target_link_libraries(nudge nudge_imgui OpenGL::GL)
	endif()

	if(MSVC)
		target_compile_options(nudge PRIVATE /std:clatest /experimental:c11atomics /Ob3)
	endif()

	# Remote debug viewer is Win32-only (uses ReadProcessMemory).
	if(WIN32)
		add_executable(nudge_viewer tools/viewer.c)
		target_link_libraries(nudge_viewer ws2_32)
	endif()

	# Copy SDL3.dll next to executables so F5 works in MSVC
	if(WIN32)
		add_custom_command(TARGET nudge POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:SDL3::SDL3-shared> $<TARGET_FILE_DIR:nudge>)
		add_custom_command(TARGET nudge_tests POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:SDL3::SDL3-shared> $<TARGET_FILE_DIR:nudge_tests>)
	endif()
endif()
