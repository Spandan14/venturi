cmake_minimum_required(VERSION 3.25)
project(Venturi LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)
add_subdirectory(tests)

target_include_directories(Venturi
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)


find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(Venturi PRIVATE Eigen3::Eigen)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(Venturi PRIVATE -Wall -Wextra -Wpedantic)
endif()
