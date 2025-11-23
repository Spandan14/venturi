add_library(engine
    sim2d.cc
    mac/mac2d.cc
    mac/mac3d.cc
)

target_include_directories(engine
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(engine PUBLIC solver utils Eigen3::Eigen)
