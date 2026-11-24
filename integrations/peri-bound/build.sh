#!/bin/sh
# Build the small Peri-side model worker against the audited peri-model source.
set -eu
src=${1:?usage: build.sh PERI_CHECKOUT EMPTY_BUILD_DIR}
out=${2:?empty build directory required}
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
mkdir -p "$out"
out=$(cd "$out" && pwd)
test -z "$(ls -A "$out")"
rev=c7d154c9985cbe08ba25e964cb7468d15014e6d0
git -C "$src" archive -o "$out/source.tar" "$rev"
tar xf "$out/source.tar" -C "$out"
mkdir -p "$out/peri-model/examples"
cp "$here/helixweave-bound.rs" "$out/peri-model/examples/"
cat "$here/single-request.rs" >> "$out/peri-model/src/anthropic/mod.rs"
(cd "$out" && bash scripts/cargo-rmcp-patched.sh build --locked --offline -p peri-model --example helixweave-bound)
