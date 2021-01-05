#!/bin/bash

curDir="$PWD"
codeDir="$curDir/src"
buildDir="$curDir/gebouw"

compiler=clang++

opts="-O3 -g -ggdb -Wall -Werror -pedantic -Wno-gnu-zero-variadic-macro-arguments -Wno-gnu-anonymous-struct -Wno-nested-anon-types -Wno-missing-braces -Wno-unused-function -Wno-writable-strings -Wno-c99-extensions -Wno-vla-extension"

mkdir -p "$buildDir"

echo "Building resampler..."

cd "$buildDir" > /dev/null

    $compiler $opts "$codeDir"/resampler.cpp -o resampler

cd "$curDir" > /dev/null
