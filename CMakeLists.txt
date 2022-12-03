# CMakeList.txt : CMake project for jabcodecpp, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("jabcodecpp")

# Add source to this project's executable.
add_library(jabcodecpp "JabCPP.cpp" "JabCPP.h" 
"src/transform.h" "src/transform.cpp" "src/sample.h" "src/sample.cpp" "src/reportError.h" "src/reportError.cpp" "src/pseudo_random.h"
"src/mask.h" "src/mask.cpp" "src/ldpc.h" "src/ldpc.cpp" "src/jabTypes.h" "src/interleave.h" "src/interleave.cpp" "src/image.h" "src/image.cpp"
"src/encoderTypes.h" "src/encoder.h" "src/encoder.cpp" "src/detectorTypes.h" "src/detector.h" "src/detector.cpp" 
"src/decoderTypes.h" "src/decoder.h" "src/decoder.cpp" "src/binarizer.h" "src/binarizer.cpp"
)

find_package(PNG REQUIRED)
target_link_libraries(jabcodecpp PRIVATE PNG::PNG)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET jabcodecpp PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
