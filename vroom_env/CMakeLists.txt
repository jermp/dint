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

add_executable(check_encoded_data check_encoded_data.cpp)
target_link_libraries(check_encoded_data
    ${Boost_LIBRARIES}
    FastPFor_lib
    streamvbyte
    MaskedVByte
    )