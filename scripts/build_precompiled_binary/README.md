This dir contains a script to build a precompiled binary for `make_prg`.

Just run from the project root:
```
scripts/build_precompiled_binary/build_precompiled_binary.sh
```
and you will find the precompiled binary at `precompiled_binary/dist/make_prg_{version}`.

It needs the container from `Dockerfile` built and properly tagged.
