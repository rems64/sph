cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17) # 14 required for make_unique; 17 required for optional
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable all warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-Wall -Wextra -Wpedantic)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(/W4)
endif()

project(sph)

file(GLOB SRC
     "src/*.cpp"
)

# ImGui source
file(GLOB IMGUI_SRC
    "dep/imgui/*.cpp"
    "dep/imgui/backends/imgui_impl_glfw.cpp"
    "dep/imgui/backends/imgui_impl_opengl3.cpp"
)

add_executable(${PROJECT_NAME} ${SRC} ${IMGUI_SRC})


add_subdirectory(dep/glad)
add_subdirectory(dep/glfw)
add_subdirectory(dep/glm)
target_link_libraries(${PROJECT_NAME} PRIVATE glad)
target_link_libraries(${PROJECT_NAME} PRIVATE glfw)
target_link_libraries(${PROJECT_NAME} PRIVATE glm::glm-header-only)
target_include_directories(${PROJECT_NAME} PRIVATE dep/glm/../)
target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_DL_LIBS})

target_include_directories(${PROJECT_NAME} PRIVATE dep/imgui/)

add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_CURRENT_SOURCE_DIR})
