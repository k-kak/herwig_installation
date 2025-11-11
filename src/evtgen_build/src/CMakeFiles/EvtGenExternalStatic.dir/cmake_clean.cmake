file(REMOVE_RECURSE
  "../lib/libEvtGenExternal.a"
  "../lib/libEvtGenExternal.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/EvtGenExternalStatic.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
