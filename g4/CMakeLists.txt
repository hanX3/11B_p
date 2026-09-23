# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
set (CMAKE_CXX_STANDARD 17)
project(HB)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find ROOT
find_package(ROOT QUIET)
if(ROOT_FOUND)
  message(STATUS "ROOT found.")
  add_definitions(-DG4ANALYSIS_USE_ROOT)
  include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "ROOT not found.")
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(HB HB.cc ${sources} ${headers})
target_link_libraries(HB ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

set(HB_SCRIPTS
  macro/gui.mac
  macro/run1.mac
  macro/run2.mac
  macro/init_vis.mac
  macro/vis.mac
  )

foreach(_script ${HB_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
