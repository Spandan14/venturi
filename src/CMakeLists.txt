add_subdirectory(engine)
add_subdirectory(solver)
add_subdirectory(utils)
add_subdirectory(renderer)
add_subdirectory(parser)

add_executable(Venturi
    main.cc
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(Venturi 
    PRIVATE 
        engine
        solver
        utils
        renderer
        parser
        Eigen3::Eigen
)
