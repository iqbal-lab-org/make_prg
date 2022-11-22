#!/usr/bin/env bash
cp ../../requirements.txt .
docker build . -t leandroishilima/make_prg_precompiled_binary_builder:0.4.0
rm -v requirements.txt
