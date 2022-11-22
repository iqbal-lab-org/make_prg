This dir contains some scripts to build a precompiled binary for `make_prg`.

You generally just need to run from the project root:
```
scripts/build_precompiled_binary/build_precompiled_binary.sh
```
and you will find the precompiled binary at `precompiled_binary/dist/make_prg_{version}`.

If for some reason you need to build the container, do:
```
cd scripts/build_precompiled_binary
./build_container.sh
```
