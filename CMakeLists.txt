cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package (SeqAn REQUIRED)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

add_executable (calculate_statistics calculate_statistics.cpp)

target_link_libraries (calculate_statistics seqan3::seqan3)
target_link_libraries (calculate_statistics ${SEQAN_LIBRARIES})
