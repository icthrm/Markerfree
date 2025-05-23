
aux_source_directory(. MRC_SRCS)

MESSAGE(STATUS "building mrc")

include_directories(${PROJECT_SOURCE_DIR})

add_library (mrc STATIC ${MRC_SRCS})

target_link_libraries(mrc)

set_target_properties(mrc PROPERTIES OUTPUT_NAME "mrc")
set_target_properties(mrc PROPERTIES CLEAN_DIRECT_OUTPUT 1)