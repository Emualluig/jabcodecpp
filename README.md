## For information about JAB Code visit the reference repository here: https://github.com/jabcode/jabcode
Jabcodecpp is simply a modified version of the original source code. All credit to the original authors.

## jabcodecpp
The goal of this project is to have a C++ version of the original library.

As well as adding a way to easily integrate jabcodecpp into other programs.

## How to use

## Changes
 - Converted code to C++20
 - Removed VLAs and flexible array member
 - Fixed array runover in placeMasterMetadataPartII caused by encodeMasterMetadata, view comments for more details, search "FIX1:"
 - Optimized deinterleaveData
 - Optimized filterBinary
 - Modified pseudo_random.h to provide a way to get the previous random number
 - Fixed potential memory leaks
 - Added writing/reading PNG to/from memory
 - Removed TIFF support (due to problems compiling to wasm)

## Building
This projects uses CMake and vcpkg to build. Currently the build script is very limited.

## Performance
