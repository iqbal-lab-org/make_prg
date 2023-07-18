# To build: docker build . -t make_prg:0.4.0
# Tagged as such, it can be used in scripts/build_precompiled_binary/build_precompiled_binary.sh to build the precompiled binary
FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update \
    && apt install -y curl graphviz graphviz-dev build-essential \
    && apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# install just
ENV JUST_VERSION="1.8.0"
RUN (curl --proto '=https' --tlsv1.2 -sSf https://just.systems/install.sh \
    | bash -s -- --tag $JUST_VERSION --to /bin) \
    && just --version

# install poetry
ENV POETRY_VERSION="1.2.2" \
    POETRY_HOME="/usr/local" \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_CREATE=false \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1
RUN (curl -sSL https://install.python-poetry.org | python3 -) \
    && poetry --version

# install make_prg
WORKDIR /make_prg
COPY . /make_prg

RUN poetry run pip install -U pip \
    && poetry install --no-ansi --only main --all-extras \
    && make_prg --version

# workaround required for pyinstaller to work
RUN cp -vr /usr/bin/* /usr/sbin/

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]