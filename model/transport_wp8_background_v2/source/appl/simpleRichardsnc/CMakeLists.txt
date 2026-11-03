# SPDX-FileCopyrightInfo: Copyright 漏 DuMux-Course contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES  params.input)

dumux_add_test(NAME paper2_transport
               SOURCES main.cc)

target_compile_features(paper2_transport PRIVATE cxx_std_17)
target_compile_options(paper2_transport PRIVATE -Wall -Wextra -Wpedantic)
