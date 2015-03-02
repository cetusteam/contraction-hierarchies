cmake_minimum_required (VERSION 2.8)
project (many2many)

find_package(Boost 1.57 COMPONENTS regex REQUIRED)
    
# create bin
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

file(GLOB MANY2MANY_SRC *.h *.cpp)

add_executable(many2many ${MANY2MANY_SRC})
target_link_libraries(many2many ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})


