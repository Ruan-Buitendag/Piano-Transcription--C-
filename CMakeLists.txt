cmake_minimum_required(VERSION 3.26)
project(C_Piano_Transcription C)

set(CMAKE_C_STANDARD 11)

add_executable(C_Piano_Transcription code/source/main.c
        code/source/dictionary.c
        code/headers/dictionary.h
        code/source/stft.c
        code/headers/stft.h
        code/source/dynamicarray.c
        code/headers/dynamicarray.h
        code/source/spectrogram.c
        code/headers/spectrogram.h
        code/source/activations.c
        code/headers/activations.h
        code/source/wav.c
        code/headers/wav.h
)


find_package(HDF5 COMPONENTS C HL)

include_directories(${HDF5_INCLUDE_DIRS} code/headers)
target_link_libraries(
        C_Piano_Transcription
        ${HDF5_LIBRARIES} # This should list all libraries.
)
