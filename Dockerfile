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

# install mafft
ARG MAFFT_VERSION="7.505"
RUN curl --proto '=https' -sSf "https://mafft.cbrc.jp/alignment/software/mafft_${MAFFT_VERSION}-1_amd64.deb" -o mafft.deb \
    && dpkg -i mafft.deb \
    && rm mafft.deb \
    && mafft --version

# Copy only requirements, to cache them in docker layer
WORKDIR /make_prg
COPY . /make_prg

RUN poetry run pip install -U pip \
  && poetry install --no-ansi --only main --all-extras

# we copy at the end so that we don't reinstall dependencies everytime there is a change
# in the code

RUN make_prg --version

SHELL ["/bin/bash", "-eo", "pipefail", "-c"]