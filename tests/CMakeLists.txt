add_executable(mac_tests_2d
    mac/mac_2d_test.cc
)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(mac_tests_2d PRIVATE engine Eigen3::Eigen)

# add_executable(mac_tests_3d
#     mac/mac_3d_test.cc
# )
#
# find_package(Eigen3 REQUIRED NO_MODULE)
# target_link_libraries(mac_tests_3d PRIVATE engine Eigen3::Eigen)

enable_testing()
add_test(NAME mac_2d_test COMMAND mac_tests_2d)
# add_test(NAME mac_3d_test COMMAND mac_tests_3d)

