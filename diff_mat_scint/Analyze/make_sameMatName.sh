#!/bin/bash

Dir="$1"
for f in $Dir/*.root; do mv "$f" "${f//_test/}"; done
