cmake_minimum_required(VERSION 2.8)
project(DS2I)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

if (UNIX)

   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

   if (USE_SANITIZERS)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
   endif ()

endif()

find_package(Boost COMPONENTS iostreams unit_test_framework filesystem system log log_setup date_time chrono REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${DS2I_SOURCE_DIR})

add_subdirectory(succinct EXCLUDE_FROM_ALL)
add_subdirectory(FastPFor EXCLUDE_FROM_ALL)

# Add streamvbyte
include_directories(streamvbyte/include)
add_library(streamvbyte STATIC streamvbyte/src/streamvbyte.c
                               streamvbyte/src/streamvbytedelta.c
)

# Add maskedvbyte
include_directories(MaskedVByte/include)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
add_library(MaskedVByte STATIC MaskedVByte/src/varintdecode.c
                               MaskedVByte/src/varintencode.c
)

add_executable(create_freq_index create_freq_index.cpp)
target_link_libraries(create_freq_index
  ${Boost_LIBRARIES}
    FastPFor_lib
    streamvbyte
    MaskedVByte
  )

add_executable(encode encode.cpp)
target_link_libraries(encode
    ${Boost_LIBRARIES}
    FastPFor_lib
    streamvbyte
    MaskedVByte
    )

add_executable(decode decode.cpp)
target_link_libraries(decode
    ${Boost_LIBRARIES}
    FastPFor_lib
    streamvbyte
    MaskedVByte
    )

#add_executable(check_encoded_data check_encoded_data.cpp)
#target_link_libraries(check_encoded_data
#    ${Boost_LIBRARIES}
#    FastPFor_lib
#    streamvbyte
#    MaskedVByte
#    )

add_executable(create_wand_data create_wand_data.cpp)
target_link_libraries(create_wand_data
  ${Boost_LIBRARIES}
  )

add_executable(queries queries.cpp)
target_link_libraries(queries
  ${Boost_LIBRARIES}
  FastPFor_lib
  streamvbyte
  MaskedVByte
  )

add_executable(dict_perf_test dict_perf_test.cpp)
target_link_libraries(dict_perf_test
  ${Boost_LIBRARIES}
  )

enable_testing()
add_subdirectory(test)
