# SPDX-FileCopyrightInfo: Copyright © DuMux-Course contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES  params.input)

dumux_add_test(NAME test_model_salt
               SOURCES main.cc)
