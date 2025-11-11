file(REMOVE_RECURSE
  "../lib/libEvtGen.a"
  "../lib/libEvtGen.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/EvtGenStatic.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
