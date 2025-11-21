add_library(engine
    engine/mac/mac2d.cc
    engine/mac/mac3d.cc
)

target_include_directories(engine
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/engine
)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(engine PRIVATE Eigen3::Eigen)

target_link_libraries(Venturi PRIVATE engine Eigen3::Eigen)
