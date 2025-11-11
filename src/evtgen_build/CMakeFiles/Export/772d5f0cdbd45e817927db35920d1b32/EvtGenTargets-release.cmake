#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "EvtGen::EvtGen" for configuration "Release"
set_property(TARGET EvtGen::EvtGen APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EvtGen::EvtGen PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libEvtGen.2.2.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libEvtGen.2.dylib"
  )

list(APPEND _cmake_import_check_targets EvtGen::EvtGen )
list(APPEND _cmake_import_check_files_for_EvtGen::EvtGen "${_IMPORT_PREFIX}/lib/libEvtGen.2.2.0.dylib" )

# Import target "EvtGen::EvtGenStatic" for configuration "Release"
set_property(TARGET EvtGen::EvtGenStatic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EvtGen::EvtGenStatic PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libEvtGen.a"
  )

list(APPEND _cmake_import_check_targets EvtGen::EvtGenStatic )
list(APPEND _cmake_import_check_files_for_EvtGen::EvtGenStatic "${_IMPORT_PREFIX}/lib/libEvtGen.a" )

# Import target "EvtGen::EvtGenExternal" for configuration "Release"
set_property(TARGET EvtGen::EvtGenExternal APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EvtGen::EvtGenExternal PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libEvtGenExternal.2.2.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libEvtGenExternal.2.dylib"
  )

list(APPEND _cmake_import_check_targets EvtGen::EvtGenExternal )
list(APPEND _cmake_import_check_files_for_EvtGen::EvtGenExternal "${_IMPORT_PREFIX}/lib/libEvtGenExternal.2.2.0.dylib" )

# Import target "EvtGen::EvtGenExternalStatic" for configuration "Release"
set_property(TARGET EvtGen::EvtGenExternalStatic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EvtGen::EvtGenExternalStatic PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libEvtGenExternal.a"
  )

list(APPEND _cmake_import_check_targets EvtGen::EvtGenExternalStatic )
list(APPEND _cmake_import_check_files_for_EvtGen::EvtGenExternalStatic "${_IMPORT_PREFIX}/lib/libEvtGenExternal.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
