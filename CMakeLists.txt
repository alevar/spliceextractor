cmake_minimum_required(VERSION 3.14)
project(spliceextractor)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

# thanks for the macro: David Grayson
macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
        endif ()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)

use_cxx11()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer -g")

include("cmake/htslib.cmake")

add_executable(spliceextractor
               src/spliceextractor.cpp
               src/arg_parse.h
               src/arg_parse.cpp)
add_dependencies(spliceextractor htslib)
target_link_libraries(spliceextractor ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)