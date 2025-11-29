add_library(fluxlang 
    parser_2d.cc)

target_include_directories(fluxlang
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${PROJECT_SOURCE_DIR}/extern/peglib
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(fluxlang PUBLIC Eigen3::Eigen)
