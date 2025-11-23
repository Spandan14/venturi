add_library(solver 
    linear_ode_solver.cc)

target_include_directories(solver
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(solver PUBLIC Eigen3::Eigen)
