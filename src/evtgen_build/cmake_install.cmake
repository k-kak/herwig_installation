# Install script for directory: /Users/user/phd/herwig_gh/src/EvtGen/R02-02-00

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/user/phd/herwig_gh")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/EvtGen")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/EvtGenBase")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/EvtGenExternal")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/EvtGenModels")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/EvtGen" TYPE FILE FILES
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/DECAY.DEC"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/DECAY.XML"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/evt.pdl"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/DECAY_2010.XML"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/DECAY_2010.DEC"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/DECAY_2009.XML"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/DECAY_2009.DEC"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/EvtGen" TYPE FILE FILES
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/AUTHORS"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/COPYING"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/EvtGen" TYPE FILE FILES
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/README.md"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/Pythia8_README.md"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/Tauola_README.md"
    "/Users/user/phd/herwig_gh/src/EvtGen/R02-02-00/History.md"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake" TYPE FILE FILES
    "/Users/user/phd/herwig_gh/src/evtgen_build/EvtGenConfig.cmake"
    "/Users/user/phd/herwig_gh/src/evtgen_build/EvtGenConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake/EvtGenTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake/EvtGenTargets.cmake"
         "/Users/user/phd/herwig_gh/src/evtgen_build/CMakeFiles/Export/772d5f0cdbd45e817927db35920d1b32/EvtGenTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake/EvtGenTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake/EvtGenTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake" TYPE FILE FILES "/Users/user/phd/herwig_gh/src/evtgen_build/CMakeFiles/Export/772d5f0cdbd45e817927db35920d1b32/EvtGenTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/EvtGen/cmake" TYPE FILE FILES "/Users/user/phd/herwig_gh/src/evtgen_build/CMakeFiles/Export/772d5f0cdbd45e817927db35920d1b32/EvtGenTargets-release.cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/user/phd/herwig_gh/src/evtgen_build/src/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/user/phd/herwig_gh/src/evtgen_build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/user/phd/herwig_gh/src/evtgen_build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
