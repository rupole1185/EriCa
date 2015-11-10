cmake_minimum_required(VERSION 2.6)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) 

project(EriCa)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules"
  )

if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()

add_subdirectory( ${EriCa_SOURCE_DIR}/OpenFOAM)
add_subdirectory( ${EriCa_SOURCE_DIR}/MariSa)

find_package(OpenMP)
find_package(MAYA REQUIRED)

include_directories( ${EriCa_SOURCE_DIR}/MariSa)
include_directories( ${MAYA_INCLUDE_DIR})

file(GLOB ErcSOURCES EriCa/*.c EriCa/*/*.c)

add_executable( EriCa ${ErcSOURCES})

target_link_libraries( EriCa m)
target_link_libraries( EriCa MariSa)
target_link_libraries( EriCa ${MAYA_LIBRARIES})

#make install to local system
install(TARGETS EriCa DESTINATION bin)

#create package
include(InstallRequiredSystemLibraries)

include(CPack)
