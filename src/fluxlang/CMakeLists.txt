add_library(fluxlang 
    parser.cc
    ast/visitor.cc
    ast/nodes.cc
    transformer.cc
    runtime.cc
)

target_include_directories(fluxlang
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${PROJECT_SOURCE_DIR}/extern/peglib
    PRIVATE
        /opt/homebrew/include
)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GLFW3 REQUIRED)

target_link_libraries(fluxlang PUBLIC glfw
    /opt/homebrew/lib/libGLEW.dylib
    Eigen3::Eigen
    imgui
)
