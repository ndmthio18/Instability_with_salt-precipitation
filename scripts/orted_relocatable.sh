#!/bin/sh
set -eu
bin_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
bundle_dir=$(dirname -- "$bin_dir")
exec "$bundle_dir/lib/ld-linux-x86-64.so.2" --library-path "$bundle_dir/lib" "$bundle_dir/bin/orted.real" "$@"
