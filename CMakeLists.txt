cmake_minimum_required(VERSION 3.21)
project(nudge C CXX)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
target_link_libraries(nudge nudge_imgui opengl32)

add_executable(nudge_tests src/test_main.c)
if(MSVC)
	target_compile_options(nudge PRIVATE /experimental:c11atomics /Ob3)
	target_compile_options(nudge_tests PRIVATE /experimental:c11atomics /Ob3)
	target_link_options(nudge_tests PRIVATE /STACK:8388608)
endif()

# Copy SDL3.dll next to executables so F5 works in MSVC
add_custom_command(TARGET nudge POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:SDL3::SDL3-shared> $<TARGET_FILE_DIR:nudge>)
add_custom_command(TARGET nudge_tests POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:SDL3::SDL3-shared> $<TARGET_FILE_DIR:nudge_tests>)
