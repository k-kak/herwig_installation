file(REMOVE_RECURSE
  "../lib/libEvtGen.2.2.0.dylib"
  "../lib/libEvtGen.2.dylib"
  "../lib/libEvtGen.dylib"
  "../lib/libEvtGen.pdb"
  ".2"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/EvtGen.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
