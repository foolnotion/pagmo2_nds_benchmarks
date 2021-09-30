cmake_minimum_required(VERSION 3.8)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# removes the need to specify -DCMAKE_EXPORT_COMPILE_COMMANDS=ON in the cli
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

project(pagmo-nds VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(./pagmo2 ./build EXCLUDE_FROM_ALL)

add_executable(nsga2_dtlz
    main.cpp
)
target_link_libraries(nsga2_dtlz PRIVATE pagmo tbb)
