dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME p1_compositional_2d
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=P1Compositional2D)
target_compile_features(p1_compositional_2d PRIVATE cxx_std_17)
target_compile_options(p1_compositional_2d PRIVATE -Wall -Wextra -Wpedantic)

dumux_add_test(NAME p1_compositional_1d
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=P1Compositional1D)
target_compile_features(p1_compositional_1d PRIVATE cxx_std_17)
target_compile_options(p1_compositional_1d PRIVATE -Wall -Wextra -Wpedantic)
