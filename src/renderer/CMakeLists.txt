add_library(renderer 
    2d/renderer_2d.cc
    3d/renderer_3d.cc
)

target_include_directories(renderer 
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
    PRIVATE
        /opt/homebrew/include
)

target_compile_definitions(renderer PRIVATE SHADERS_3D="${CMAKE_CURRENT_SOURCE_DIR}/3d/shaders/")

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenGL REQUIRED)
find_package(GLFW3 REQUIRED)
# find_package(GLEW REQUIRED)

# find the slang library
find_library(SLANG_LIB
    NAMES slang
    PATHS /opt/homebrew/lib
    REQUIRED
)

target_include_directories(renderer PRIVATE ${SLANG_INCLUDE_DIR})

target_link_libraries(renderer PUBLIC 
    ${SLANG_LIB}
    OpenGL::GL
    glfw
    /opt/homebrew/lib/libGLEW.dylib
    Eigen3::Eigen
    imgui
)
