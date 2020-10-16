cmake_minimum_required(VERSION 3.3)
project(E2R)



find_package(Geant4 REQUIRED)

if(NOT Geant4_Found)
	message("Geant4 is not found.")
endif()

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
include_directories(${PROJECT_SOURCE_DIR}/include)


add_executable(main.cc ${sources} ${headers})


