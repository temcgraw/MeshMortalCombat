cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


project(DestructibleVoxelSystem)

# Compile dependencies
add_subdirectory(./3rd_party/glfw-3.3.8)
add_subdirectory(./3rd_party/assimp)

set(IMGUI_FILES
./3rd_party/imgui/imgui.cpp
./3rd_party/imgui/imgui_draw.cpp
./3rd_party/imgui/imgui_tables.cpp
./3rd_party/imgui/imgui_widgets.cpp
./3rd_party/imgui/imgui_demo.cpp
./3rd_party/imgui/backends/imgui_impl_glfw.cpp # platform
./3rd_party/imgui/backends/imgui_impl_opengl3.cpp # render
)

include_directories(./
                    ./3rd_party
                    ./3rd_party/glad/include
                    ./3rd_party/imgui
                    ./3rd_party/imgui/backends
		    ./3rd_party/stb_image)


file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS source/*.cpp source/*.h)

list(APPEND SOURCES
    3rd_party/glad/src/glad.c  # glad files
)

add_executable(${PROJECT_NAME} ${SOURCES} ${IMGUI_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE include)

target_link_libraries(${PROJECT_NAME} PUBLIC glfw assimp)



# set output directories
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE   ${CMAKE_CURRENT_SOURCE_DIR}/bin/Release)
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG     ${CMAKE_CURRENT_SOURCE_DIR}/bin/Debug)
if(WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin/Debug)
elseif(UNIX)
endif()

# Copy resources and shaders to target directory
# set the current build output directory
set(OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin/$<CONFIG>)


add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/resource
    ${OUTPUT_DIR}/resource
    COMMAND ${CMAKE_COMMAND} -E echo
    "copying ${CMAKE_CURRENT_SOURCE_DIR}/resource to ${OUTPUT_DIR}/resource"
    COMMENT "Copying resources to build directory")

add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/shaders
    ${OUTPUT_DIR}/shaders
    COMMAND ${CMAKE_COMMAND} -E echo
    "copying ${CMAKE_CURRENT_SOURCE_DIR}/shaders to ${OUTPUT_DIR}/shaders"
    COMMENT "Copying shaders to build directory")
add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/outputs
    ${OUTPUT_DIR}/outputs
    COMMAND ${CMAKE_COMMAND} -E echo
    "copying ${CMAKE_CURRENT_SOURCE_DIR}/outputs to ${OUTPUT_DIR}/outputs"
    COMMENT "Copying output folder to build directory")
