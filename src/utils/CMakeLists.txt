add_library(utils INTERFACE)

target_include_directories(utils
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(utils INTERFACE Eigen3::Eigen)
