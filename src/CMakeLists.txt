add_executable(create_freq_index create_freq_index.cpp)
target_link_libraries(create_freq_index
  ${Boost_LIBRARIES}
    FastPFor_lib
    streamvbyte
    MaskedVByte
  )

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

add_executable(pair_wise_intersect pair_wise_intersect.cpp)
target_link_libraries(pair_wise_intersect
  ${Boost_LIBRARIES}
  FastPFor_lib
  streamvbyte
  MaskedVByte
  )

add_executable(check_index check_index.cpp)
target_link_libraries(check_index
  ${Boost_LIBRARIES}
  FastPFor_lib
  streamvbyte
  MaskedVByte
  )

add_executable(dict_perf_test dict_perf_test.cpp)
target_link_libraries(dict_perf_test
  ${Boost_LIBRARIES}
  )