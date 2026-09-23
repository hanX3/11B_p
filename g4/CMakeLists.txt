# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(HB)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default.
# Set WITH_GEANT4_UIVIS=OFF to build a batch-mode-only executable.
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions.
include(${Geant4_USE_FILE})

# Find ROOT.
find_package(ROOT QUIET)
if(ROOT_FOUND)
  message(STATUS "ROOT found.")
  add_definitions(-DG4ANALYSIS_USE_ROOT)
else()
  message(FATAL_ERROR "ROOT not found.")
endif()

# Project include directories.
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${Geant4_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
)

# Source files.
file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

# Executable.
add_executable(HB HB.cc ${sources} ${headers})
target_link_libraries(HB ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Runtime files

# Keep batch.sh in the build directory.
configure_file(
  ${PROJECT_SOURCE_DIR}/batch.sh
  ${PROJECT_BINARY_DIR}/batch.sh
  COPYONLY
)

# Copy the complete macros directory to the build directory at build time.
# This ensures newly added validation macros are also available under build/macros/.
add_custom_target(copy_macros ALL
  COMMAND ${CMAKE_COMMAND} -E remove_directory
          ${PROJECT_BINARY_DIR}/macros
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${PROJECT_SOURCE_DIR}/macros
          ${PROJECT_BINARY_DIR}/macros
  COMMENT "Copying macros directory to build directory"
)

# Ensure macros are copied when building the main executable.
add_dependencies(HB copy_macros)
