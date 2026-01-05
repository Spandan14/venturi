add_library(imgui 
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/backends/imgui_impl_opengl3.cpp
    imgui/backends/imgui_impl_glfw.cpp
    imgui/implot.cpp
    imgui/implot_items.cpp
)

target_include_directories(imgui PUBLIC 
    imgui
    imgui/backends
)

find_package(GLFW3 REQUIRED)

target_link_libraries(imgui PUBLIC glfw)
