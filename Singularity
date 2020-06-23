Bootstrap: docker
from: python:3.8-alpine

%environment
    export PATH=/usr/local:/usr/local/bin:$PATH

%post
    cd / || exit 1

    # INSTALL DEPS
    apk update && apk upgrade && \
    apk --no-cache add lapack libstdc++ openblas-dev bash fd && \
    apk --no-cache add --virtual .builddeps g++ gcc gfortran musl-dev lapack-dev git

    # CLONE
    URL="https://github.com/rmcolq/make_prg.git"
    DIR="make_prg"
    TAG="0bb4a27af50b9d8cd80995fca9d93ed4b8f78f50"
    git clone "$URL" "/${DIR}"
    cd "/${DIR}" || exit 1
    git checkout "$TAG"

    # INSTALL
    python -m pip install --no-cache-dir .

    # remove bulky cache
    cd / || exit 1
    apk del .builddeps && rm -rf /root/.cache

%test
    export PATH=/usr/local:/usr/local/bin:$PATH
    cd /make_prg
    python -m pip install --no-cache-dir nose hypothesis
    nosetests tests/
    make_prg --help


%labels
    Author Michael Hall
    Version 0bb4a27
    Website https://github.com/rmcolq/make_prg


%help
    This is a container for running the make_prg python command-line tool only.

    # Example usage
    singularity exec <container> make_prg --help
