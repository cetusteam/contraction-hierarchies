cmake_minimum_required (VERSION 3.1)
project (contraction-hierarchies)

find_package(Boost 1.57 COMPONENTS regex REQUIRED)
find_package(Protobuf)

add_subdirectory(protos/ch)

# create bin
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)
set (GENERATED_DIR ${CMAKE_SOURCE_DIR}/generated)

file(GLOB_RECURSE INCLUDE_FILES
    stats/*.h
    command/*.h
    io/*.h
    datastr/*.h
    processing/*.h
)

file(GLOB CHMAIN_FILES *.h *.cpp)

include_directories(
	${PROTOBUF_INCLUDE_DIRS}
	${BOOST_INCLUDE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}
)

set(CHMAIN_SRC
    ${CHMAIN_FILES}
    ${INCLUDE_FILES}
)

# This program is capable to generate a graph required to run many2many
add_executable(node-order ${CHMAIN_SRC})
# We don't use the "-e" flag yet
# set_target_properties(node-order PROPERTIES COMPILE_FLAGS "-DCOUNT_SHORTCUT_ORIGINAL_EDGES")
target_link_libraries(node-order
    protos-ch
    ${Boost_LIBRARIES}
    ${PROTOBUF_LIBRARIES}
)

# Disabled for now because we don't generate paths, yet.
# add_executable(contraction ${CHMAIN_SRC})
# For many2many we don't need this parameter. However, it's required to
# calculate the paths
# set_target_properties(contraction PROPERTIES COMPILE_FLAGS "-DUSE_CH_EXPAND")
# target_link_libraries(contraction ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})

add_subdirectory(many)
