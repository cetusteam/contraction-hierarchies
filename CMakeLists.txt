cmake_minimum_required (VERSION 2.8)
project (contraction-hierarchies)

find_package(Boost 1.57 COMPONENTS regex REQUIRED)
    
# create bin
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

file(GLOB PROTO_FILES ${SOLVER_PROTOS}/*.proto)

FOREACH(PROTOFILE ${PROTO_FILES})
    GET_FILENAME_COMPONENT(FILE_WE ${PROTOFILE} NAME_WE)
    LIST(APPEND GENERATED_PROTO_FILES "${GENERATED_DIR}/${FILE_WE}.pb.cc")
    LIST(APPEND GENERATED_PROTO_FILES "${GENERATED_DIR}/${FILE_WE}.pb.h")
ENDFOREACH()

FOREACH(GENERATED_PROTO_FILE ${GENERATED_PROTO_FILES})
    set_source_files_properties(${GENERATED_PROTO_FILE} PROPERTIES GENERATED TRUE )
ENDFOREACH()


function(ADD_PROTO_FILES_VARS)
    ADD_CUSTOM_COMMAND(
       OUTPUT ${GENERATED_PROTO_FILES}
       COMMAND mkdir -p ${CMAKE_SOURCE_DIR}/generated
       COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} --cpp_out=${GENERATED_DIR} --proto_path=${SOLVER_PROTOS} ${SOLVER_PROTOS}/*.proto
       DEPENDS ${PROTO_FILES}
       )
endfunction()


file(GLOB_RECURSE SUBDIR_INCLUDE_FILES
    stats/*.h
    command/*.h
    io/*.h
    datastr/*.h
    processing/*.h
)

file(GLOB CHMAIN_FILES *.h *.cpp)

set(CHMAIN_SRC
    ${CHMAIN_FILES}
    ${SUBDIR_INCLUDE_FILES}
)


# This program is capable to generate a graph required to run many2many
add_executable(node-order ${CHMAIN_SRC})
# We don't use the "-e" flag yet
# set_target_properties(node-order PROPERTIES COMPILE_FLAGS "-DCOUNT_SHORTCUT_ORIGINAL_EDGES")
target_link_libraries(node-order ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})

# Disabled for now because we don't generate paths, yet.
# add_executable(contraction ${CHMAIN_SRC})
# For many2many we don't need this parameter. However, it's required to
# calculate the paths
# set_target_properties(contraction PROPERTIES COMPILE_FLAGS "-DUSE_CH_EXPAND")
# target_link_libraries(contraction ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})


add_subdirectory(many)
