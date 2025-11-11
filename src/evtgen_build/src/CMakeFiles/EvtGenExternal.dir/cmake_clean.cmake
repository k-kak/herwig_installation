file(REMOVE_RECURSE
  "../lib/libEvtGenExternal.2.2.0.dylib"
  "../lib/libEvtGenExternal.2.dylib"
  "../lib/libEvtGenExternal.dylib"
  "../lib/libEvtGenExternal.pdb"
  ".2"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/EvtGenExternal.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
